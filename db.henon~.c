/* 
	AUTHOR:			Daniel Bennett
	DATE:			10/01/2015
	DESCRIPTION:	Framework for recursive oscillators
					First Attempt - using henon
		
	Version:		1.0
	
	fixes this version
		-- 64bit MAX conversion
		-- Adapt to accept float input
		-- Remove debugging and make Bang output the current values in the array to the right outlet 
			(so we can copy and replay sections by forcing initial conditions!)
				(freezing and replaying chaos foregrounds the determinism in noise - play this against prepared speaker)
		
		-- Simplified calculations - time coordinates for outputs not required
		-- tidied code
	
		-- Scaling from 0-127 to useable ranges for variables (allows standardisation across range)
		-- Add limits on incoming numbers	
		-- Prevent numbers heading to infinity 
		-- Add method to change interpolation type on receiving symbol "interpolation" followed by an int
		-- Add method to turn on/off clipping of output to range -1 ... 1
		-- Add method to re-initialise array with random values on receiving "reset"

	Issues:
		--	While iteration rate is variable, I treat it as constant for the sake of interpolation,
			indexing current position between outputs at n-2 and n-1 with phase. This is obviously 
			going to result in distortion when pitch shifts, but since there are no real values between 
			integer values of n to deviate from; we cannot guess the changes in iteration rate that are 
			coming before the next iteration, and the signal is anyway noisy in most areas, I'm fine with this.

	TODO:
		Repeat for other algorithms

	ASPIRATIONAL TODO:	
		Change interpolation type handling to use function pointer
		Better workaround for preventing numbers jumping to infinity?
		Find a way of moving around the *interesting* parameter space - ie where no repetition borders on cycles, 
			and away from steady states/out of bounds


*/

#include "ext.h"
#include "z_dsp.h"
#include "ext_obex.h"
#include <time.h>
#include <stdlib.h>

#define OUTPUT_COUNT 4  // store enough values to calculate next value and do interpolation

// nice useable range for equation
#define MIN_A -2.0
#define MAX_A -1.0
#define MIN_B -0.25
#define MAX_B 0.25

// nice useable range for equation v2
//#define MIN_A -1.78
//#define MAX_A -0.3
//#define MIN_B -0.3
//#define MAX_B 0.88

#define SCALEOUTPUT 0.6


//Structure that defines the max object
typedef struct _henon {
	t_pxobject obj;						//"Header for any non-ui signal processing object"
	void *outlet2;						// pointer to the 2nd outlet - list out
	double outputs[OUTPUT_COUNT];		// last x output values from algorithm
	double sl;							// sample length (storing reduces no of calculations neccessary)
	double iter_rate;					// current iteration rate 
	double phase;						// current phase
	double a;							// equation variable
	double b;							// equation variable
	t_int interpolation_type;				// 0=none, 1=linear, 2-bspline
	t_int clip;							// clipping of output to -1...1, 1/0:on/off 
	t_int iter_conn;					// is signal connected to iteration inlet
	t_int a_conn;						//etc.
	t_int b_conn;						//etc.

}	t_henon;


static t_class *henon_class;	// pointer to the class of this object

// function prototypes
// MSP infrastructure functions
void	*henon_new(t_symbol *s, short argc, t_atom *argv);						// create new instance
void	henon_dsp(t_henon *x, t_signal **sp, short *count);					// connect to dsp
void	henon_dsp64(t_henon *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void	henon_assist(t_henon *x, void *b, long msg, long arg, char *dst);		// explain inputs

// handle incoming symbols
void	henon_float(t_henon *x, double f);										// handle incoming float
void	henon_bang(t_henon *x, double f);											// handle incoming bang
void	henon_reset(t_henon *x, t_symbol *msg, short argc, t_atom *argv);			// handle incoming "reset" string
void	henon_clip(t_henon *x, t_symbol *msg, short argc, t_atom *argv);			// handle incoming "clip" string
void	henon_ip_type(t_henon *x, t_symbol *msg, short argc, t_atom *argv);		// handle incoming "interpolation" string

// My infrastructure functions
double	infr_scale_param(double in, double min, double max);		// scales incoming floats from 0 - 127 to range specified

// Audio Calc functions
t_int	*henon_perform(t_int *w);												// do the bit that makes sound
void henon_perform64(t_henon *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);

double	au_henon_calc(double prev, double prevprev, double a, double b);		// calculates the next henon value
double	au_interpolate_lin(double outputs[], double x);							// interpolates value at current sample position
double	au_interpolate_bspline(double outputs[], double x);						// interpolates value at current sample position


/************************************************************

!!!!!!!!!!!!	MSP INFRASTRUCTURE FUNCTIONS	!!!!!!!!!!!!

*************************************************************/


// initialization routine 
int main (void)
{
	henon_class = class_new("db.henon~", (method)henon_new, (method)dsp_free, 
		sizeof(t_henon), 0L, A_GIMME, 0);

	// register methods to handle incoming messages
	class_addmethod(henon_class, (method)henon_dsp, "dsp", A_CANT, 0);		// Old 32-bit MSP dsp chain compilation for Max 6
	class_addmethod(henon_class, (method)henon_dsp64,	"dsp64", A_CANT, 0);		// New 64-bit MSP dsp chain compilation for Max 6
	class_addmethod(henon_class, (method)henon_assist, "assist", A_CANT, 0);
	class_addmethod(henon_class, (method)henon_float, "float", A_FLOAT, 0);
	class_addmethod(henon_class, (method)henon_bang, "bang", A_FLOAT, 0);
	class_addmethod(henon_class, (method)henon_reset, "reset", A_GIMME, 0);
	class_addmethod(henon_class, (method)henon_clip, "clip", A_GIMME, 0);
	class_addmethod(henon_class, (method)henon_ip_type, "interpolation", A_GIMME, 0);

	class_dspinit(henon_class);
	class_register(CLASS_BOX, henon_class);

	post("db.henon~ by Daniel Bennett skjolbrot@gmail.com");
	post("Sonification of the Henon Map with variable iteration rate and interpolation");
	post("args: 1: Iteration rate (float) 2: param a (float 0-127), 3:param b (float, 0-127), ");
	post("4: interpolation type (0 = none, 1 = linear, 2 = bspline(default), 5: clip signal to -1/1 (1 = on, 0 = off(default)) ");
	// report to the MAX window
	return 0;
}

// function to create new instance and initialise its parameters
void *henon_new(t_symbol *s, short argc, t_atom *argv)
{
	float iter = 8400, a = 20 , b = 20;
	long interp = 2, clip = 0;
	short i, j;

	t_henon *x = object_alloc(henon_class); // set aside memory for the struct for the object
	dsp_setup((t_pxobject *)x, 3); // call routine dsp_setup, connecting this object to the dsp chain. 3 inlets

	x->outlet2 = listout((t_object *)x); // add a list outlet 
	outlet_new((t_object *)x, "signal"); // add a signal outlet 



	atom_arg_getfloat(&iter, 0, argc, argv); // MaxMSP function, copies value from arguments on numbered max object to pointer specified
	atom_arg_getfloat(&a, 1, argc, argv); 
	atom_arg_getfloat(&b, 2, argc, argv); 
	atom_arg_getlong(&interp, 3, argc, argv); 
	atom_arg_getlong(&clip, 4, argc, argv); 

	// initialise values in object struct
	x->phase = 0;
	x->sl = 1.0 / sys_getsr();
	x->iter_rate = iter;
	x->a = a;
	x->b = b;
	x->interpolation_type = (short) interp;
	x->clip = (short) clip;

	// initialize output array
	for (i = OUTPUT_COUNT-1, j = 0; i >= 0; i--, j++) {
		x->outputs[j]  = j / 10.0;
	}
	return x;
}


//function to connect to DSP chain
void henon_dsp(t_henon *x, t_signal **sp, short *count)
{
// Check sample rate & reallocate memory if so
//	int i;
	//get status of inlets
	x->iter_conn = count[0];
	x->a_conn = count[1];
	x->b_conn = count[2];
	x->phase = 0;
	// get period of sampling
	x->sl = 1.0 / sys_getsr();
	
	// dsp add method
	dsp_add(henon_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, 
		sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}


// this is the Max 6 version of the dsp method -- it registers a function for the signal chain in Max 6,
// which operates on 64-bit audio signals.
void henon_dsp64(t_henon *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{	
	// instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
	// the arguments passed are:
	// 1: the dsp64 object passed-in by the calling function
	// 2: the symbol of the "dsp_add64" message we are sending
	// 3: a pointer to your object
	// 4: a pointer to your 64-bit perform method
	// 5: flags to alter how the signal chain handles your object -- just pass 0
	// 6: a generic pointer that you can use to pass any additional data to your perform method
	
		//get status of inlets
	x->iter_conn = count[0];
	x->a_conn = count[1];
	x->b_conn = count[2];
	x->phase = 0;
	// get period of sampling
	x->sl = 1.0 / sys_getsr();

	object_method(dsp64, gensym("dsp_add64"), x, henon_perform64, 0, NULL);
}


// assist messages for the inlets and outlets. Hovering causes MAX to send msg "ASSIST_INLET" or "ASSIST OUTLET" (or their enumeration by the look of it) to the object, 
// along with an "arg" indicating the number of the inlet/outlet
void henon_assist(t_henon *x, void *b, long msg, long arg, 
	char *dst)
{
	if (msg==ASSIST_INLET){
		switch (arg) {
		case 0: sprintf(dst,"(signal/float) iteration rate - frequency for generation of new values from equation"); break;
		case 1: sprintf(dst,"(signal/float) alpha (equation variable)"); break;
		case 2: sprintf(dst,"(signal/float) beta (equation variable)"); break;
		default: break;
		}
	}
}

/************************************************************
!!!!!!!!!!!!	INCOMING MESSAGE HANDLING		!!!!!!!!!!!!
*************************************************************/

// MSG float input, called whenever float is sent to ANY input
void henon_float(t_henon *x, double f)
{
	int inlet = ((t_pxobject*)x)->z_in;
	switch(inlet){
		case 0: x->iter_rate = f; break;
		case 1: x->a = f; break;
		case 2: x->b = f; break;
		default: break;
	}
}

// MSG BANG input - outputs list of values in output buffer
void henon_bang(t_henon *x, double f)
{
	t_atom outlist[OUTPUT_COUNT]; 
	short i; 

	for (i=0; i < OUTPUT_COUNT ; i++) { 
		atom_setfloat(outlist+i,*(x->outputs)+i);
	} 

	outlet_list(x->outlet2, 0L, OUTPUT_COUNT, outlist); 

}

// MSG "reset" symbol input - reinitialises outputs, either resets to default, or can specify list
void	henon_reset(t_henon *x, t_symbol *msg, short argc, t_atom *argv)		
{	
	short i;
	if(argc >= 1){
		for (i = 1 ; i <= OUTPUT_COUNT-1 && i < argc; i++) {
		x->outputs[i]  = atom_getfloat(argv + i);
		}
	} else {

	for (i = 0 ; i <= OUTPUT_COUNT-1 ; i++) {
		x->outputs[i]  = i / 10.0;
		}
	}
}
// MSG "clip" symbol input, turns on/off clipping to -1...1
void	henon_clip(t_henon *x, t_symbol *msg, short argc, t_atom *argv)
{
	if(argc >= 1){
		x->clip = atom_getintarg(0,argc,argv);
	}
}

// MSG "clip" symbol input, turns on/off clipping to -1...1
void	henon_ip_type(t_henon *x, t_symbol *msg, short argc, t_atom *argv)
{
	if(argc >= 1){
		x->interpolation_type = atom_getintarg(0,argc,argv);
	}
}


/************************************************************

!!!!!!!!!!!!	MY HELPER FUNCTIONS		!!!!!!!!!!!!

*************************************************************/

// scales float in range 0 - 127 to float in range min - max
double infr_scale_param(double in, double min, double max)
{
	//Clip incoming values to range
	if(in < 0.0) { in = 0.0;}
	if(in > 127.0) { in = 127.0;}

	//Scale input & return
	return min + (in * (max - min) / 127.0 );
}

/************************************************************

!!!!!!!!!!!!	AUDIO CALC FUNCTIONS		!!!!!!!!!!!!

*************************************************************/

// linear interpolation 
// returns pos between values n-2 and n-1 (defined as 0...1 on x axis) indexed by current phase
double au_interpolate_lin(double outputs[], double x)
{
	// start out just doing linear interpolation
	int n = OUTPUT_COUNT - 1;
	return outputs[n-1]  + x*(outputs[n-2] - outputs[n-1] );
}

// b spline interpolation 
// returns pos between values n-2 and n-1 (defined as 0...1 on x axis) indexed by current phase
double au_interpolate_bspline (double outputs[], double x)
{
	int n = OUTPUT_COUNT - 1;

	// 4-point, 3rd-order B-spline (x-form) (lightly adapted from Olli Niemitalo's "pink elephant" paper)
	double ym1py1 = outputs[n-3] + outputs[n-1] ;
	double c0 = 1/6.0*ym1py1 + 2/3.0*outputs[n-2] ;
	double c1 = 1/2.0*(outputs[n-1] - outputs[n-3] );
	double c2 = 1/2.0*ym1py1 - outputs[n-2] ;
	double c3 = 1/2.0*(outputs[n-2] - outputs[n-1] ) + 1 / 6.0 * (outputs[n] - outputs[n-3] );
	return ((c3*x+c2)*x+c1)*x+c0;
}


// henon map calculation y[n] = a * y[n-1]^2 + b * y[n-2] + 1
double au_henon_calc(double prev, double prevprev, double a, double b)
{
	double next;
	next =  (a * pow(prev,2)) + (b * prevprev) + 1,1000;

	// protect against areas of map where output shoots off to infinity
	if(next > 100 || next < -100) { next = 0.0;}

	return next;
}



/************************************************************

!!!!!!!!!!!!		MAX PERFORM ROUTINE			!!!!!!!!!!!!

*************************************************************/

t_int *henon_perform(t_int *w)
{	// dereference the object and declare variables to hold contents
	// corresponding to the structure declared in dsp_add above, in the vector received from the dsp chain 
	t_henon *x = (t_henon *) (w[1]);	//1 is pointer to object
	t_float *iter = (t_float *) (w[2]);			//2 is in1 (iter rate)
	t_float *a = (t_float *) (w[3]);			//3 is in2	(a)
	t_float *b = (t_float *) (w[4]);			//4 is in3	(b)
	t_float *waveout = (t_float *) (w[5]);		//5 is the wave output 
	t_int sample = w[6];						//6 is the number of samples per vector

	// variables for handling float instead of signal inputs 
	t_int iter_conn = x->iter_conn;
	t_int a_conn = x->a_conn;
	t_int b_conn = x->b_conn;
	double s_iter, s_a, s_b;		// value to hold connected inlet values this sample
	double out = 0.0 ;

	// other variables
	double holdint;
	int n = OUTPUT_COUNT - 1;
	int i;
	
	// check for mute~ - if disabled, jump immediately to the next object in the dsp chain
	if (x->obj.z_disabled) 
		return w + 7;

	//CALCS *****************************************************************************************
	while(sample--){ // for each location in the audio vector

		// Pick signal values if selected, otherwise float
		s_iter = iter_conn ? *iter++ : x->iter_rate;
		s_a = a_conn ? *a++ : x->a;
		s_b = b_conn ? *b++ : x->b;

		// scale inputs
		s_a = infr_scale_param(s_a, MIN_A, MAX_A);
		s_b = infr_scale_param(s_b, MIN_B, MAX_B); 

		// calculate new phase 
		x->phase = x->phase + fabs(s_iter * x->sl);	//signal input on iteration rate

		// If phase wraps increment array queue and calculate new value 
		if( x->phase >= 1.0 ){

			// set phase to non-integer portion of phase
			x->phase = modf(x->phase, &holdint); 

			//  move all values in array back 1
			for (i = 0; i <= n-1; i++){
				x->outputs[i]  = x->outputs[i+1] ;
			}
			// calculate new henon value
				x->outputs[n]  =  au_henon_calc(x->outputs[n-1] , x->outputs[n-2] , s_a, s_b);
		} 


		// UPDATE OUTPUTS
		// do interpolation to get value for this sample, then scale to useable range

		switch (x->interpolation_type) {
			case 0: //no interpolation, just take last value in array
				out = x->outputs[OUTPUT_COUNT - 1]  * SCALEOUTPUT;
				break;
			case 1: //linear interpolation, 
				out = au_interpolate_lin(x->outputs, x->phase) * SCALEOUTPUT;
				break;
			case 2: //bspline interpolation
				out = au_interpolate_bspline(x->outputs, x->phase) * SCALEOUTPUT; // bspline is quieter, compensate
				break;
			default: // catch case - linear interpolation (should never happen ... but!)
				out = au_interpolate_lin(x->outputs, x->phase) * SCALEOUTPUT;
				break;
		}

		// if clip != 0, do clipping
		if (x->clip) {
		out = out < -1.0 ?  -1.0 : out;
		out = out > 1.0 ?  1.0 : out;
		}

		// output
		*waveout++ = (t_float)out;

	}

	// CLOSE BUSINESS FOR THIS VECTOR
	return w + 7;		// 6 items in our dsp vector, so 8th will be the next object in the chain. Return pointer to it.	
}





// this is 64-bit perform method for Max 6
void henon_perform64(t_henon *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
// dereference the object and declare variables to hold contents
	// corresponding to the structure declared in dsp_add above, in the vector received from the dsp chain 
	t_double *iter = (t_double *) (ins[0]);			//2 is in1 (iter rate)
	t_double *a = (t_double *) (ins[1]);			//3 is in2	(a)
	t_double *b = (t_double *) (ins[2]);			//4 is in3	(b)
	t_double *waveout = (t_double *) (outs[0]);		//5 is the wave output 
	t_int sample =sampleframes;						//6 is the number of samples per vector

	// variables for handling float instead of signal inputs 
	t_int iter_conn = x->iter_conn;
	t_int a_conn = x->a_conn;
	t_int b_conn = x->b_conn;
	double s_iter, s_a, s_b;		// value to hold connected inlet values this sample
	double out = 0.0 ;

	// other variables
	double holdint;
	int n = OUTPUT_COUNT - 1;
	int i;
	
	// check for mute~ - if disabled, jump immediately to the next object in the dsp chain
	if (x->obj.z_disabled) 
		return ;

	//CALCS *****************************************************************************************
	while(sample--){ // for each location in the audio vector

		// Pick signal values if selected, otherwise float
		s_iter = iter_conn ? *iter++ : x->iter_rate;
		s_a = a_conn ? *a++ : x->a;
		s_b = b_conn ? *b++ : x->b;

		// scale inputs
		s_a = infr_scale_param(s_a, MIN_A, MAX_A);
		s_b = infr_scale_param(s_b, MIN_B, MAX_B); 

		// calculate new phase 
		x->phase = x->phase + fabs(s_iter * x->sl);	//signal input on iteration rate

		// If phase wraps increment array queue and calculate new value 
		if( x->phase >= 1.0 ){

			// set phase to non-integer portion of phase
			x->phase = modf(x->phase, &holdint); 

			//  move all values in array back 1
			for (i = 0; i <= n-1; i++){
				x->outputs[i]  = x->outputs[i+1] ;
			}
			// calculate new henon value
				x->outputs[n]  =  au_henon_calc(x->outputs[n-1] , x->outputs[n-2] , s_a, s_b);
		} 


		// UPDATE OUTPUTS
		// do interpolation to get value for this sample, then scale to useable range

		switch (x->interpolation_type) {
			case 0: //no interpolation, just take last value in array
				out = x->outputs[OUTPUT_COUNT - 1]  * SCALEOUTPUT;
				break;
			case 1: //linear interpolation, 
				out = au_interpolate_lin(x->outputs, x->phase) * SCALEOUTPUT;
				break;
			case 2: //bspline interpolation
				out = au_interpolate_bspline(x->outputs, x->phase) * SCALEOUTPUT; // bspline is quieter, compensate
				break;
			default: // catch case - linear interpolation (should never happen ... but!)
				out = au_interpolate_lin(x->outputs, x->phase) * SCALEOUTPUT;
				break;
		}

		// if clip != 0, do clipping
		if (x->clip) {
		out = out < -1.0 ?  -1.0 : out;
		out = out > 1.0 ?  1.0 : out;
		}

		// output
		*waveout++ = out;

	}
}

