#include "mex.h"
#include "spa.h"
#ifdef _WIN32
#include <string.h>
#define strcasecmp stricmp
#else
#include <strings.h>
#endif

//#define STRUCT_OUT

/*
 * siSunPosition_mex is a MEX interface to NREL's solar position algorithm, saved in spa.{h,c}.
 *
 * calling conventions:
 * [zenith, azimuth, earth_sun_dist] = siSunPosition_mex(time,position);
 *		with time a UTC datenum, and position a struct with fields latitude, longitude, and altitude (altitude in meters)
 *
 * [...] = siSunPosition_mex(time, latitude, longitude, altitude);
 *		to directly specify rather than using a structure
 *
 * [~,~,~, sunrise, sunset, suntransit] = siSunPosition_mex(...);
 *		to also calculate sunrise/sunset/sun transit (i.e. solar noon) times.  Outputs are fractional hours from the start of the day, accurate to nearest 30 seconds
 *		NOTE: this causes a more extended version of the calculation, which may slow things down slightly for a large number of points
 *
 * [...] = siSunPosition_mex(time, ..position_spec.. , timezone);
 *		lets you specify a numeric timezone in hours east of GMT (negative is west).
 *		Most of the code works in UTC, but it may be desirable to use this form when getting sunrise/sunset/sun transit times, since otherwise the SPA may get the sunrise/sunset time for the wrong day
 *
 * [...] = siSunPosition_mex(time, ..pos_spec.. , type_str);
 *		lets you specify the calculation type done.  This shouldn't really be necessary, since right now there are only two types that you can get the data out from anyway, and the best one will be selected automatically depending on how many outputs you request.
 *		valid options are the strings that correspond to the enumeration for this in spa.h:
 *			SPA_ZA, SPA_ZA_INC, SPA_ZA_RTS, SPA_ALL
 */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
/* variable declarations here */
spa_data sData;
int i,n,m;
int typeField;
mxArray *times, **times_o;
mxArray * pos = NULL;
double *t,*mon,*day,*hr,*min,*sec;
double *zenith, *azimuth, *es_dist, *sunrise, *sunset, *suntransit;
#ifdef STRUCT_OUT
mxArray *mZen, *mAz, *mR, *mT;
double *T_, *T_o;
const char *ofns[] = { "time", "zenith", "azimuth", "earthsundistance" };
#endif

/* code here */
if(nrhs < 2 || (nrhs < 4 && !mxIsStruct(prhs[1])) && !mxIsClass(prhs[1],"bu.science.geography.Position")) {
	mexErrMsgIdAndTxt("siSunPosition:input","no position specified");
}
//Get times in datevec form for passing into the solar position algorithm
times_o = (mxArray **)prhs;
mexCallMATLAB(1,&times,1,times_o,"datevec");
t = mxGetPr(times);
//Set some basics in the sData structure
sData.timezone = 0;
sData.delta_t = 32.184+35/*TAI-UTC*/+0.3/*UT1-UTC*/; //values from http://maia.usno.navy.mil/ser7/ser7.dat on Nov 1 2012
//really, should be getting TAI-UTC from http://maia.usno.navy.mil/ser7/tai-utc.dat
//UT1-UTC is always by definition < 0.9 seconds, since otherwise a leap second will be inserted/removed  For our precision purposes, we can probably just ignore this, but if we needed that kind of precision, we'd want to look it up from http://maia.usno.navy.mil/ser7/mark3.out
sData.pressure = 1000;//unless we want to be really fancy
sData.temperature = 15;
sData.slope = 0;
sData.azm_rotation = 0;
sData.atmos_refract = 0.5667; //typical value according to spa.h
sData.function = SPA_ZA; //just calculate zenith and azimuth

// Set position based on inputs
if( mxIsStruct(prhs[1]) ) {
	pos = (mxArray *)prhs[1];
} else if ( mxIsClass(prhs[1],"bu.science.geography.Position") ) {
	mexCallMATLAB(1,&pos,1,(mxArray **)(&prhs[1]),"struct");
}
if( pos ) {
	mxArray * tmpVar;
	tmpVar = mxGetField(pos, 0, "latitude");
	if(tmpVar) { sData.latitude = mxGetScalar(tmpVar); } else { mexErrMsgIdAndTxt("siSunPosition:input","no latitude specified in struct position"); }
	tmpVar = mxGetField(pos, 0, "longitude");
	if(tmpVar) { sData.longitude = mxGetScalar(tmpVar); } else { mexErrMsgIdAndTxt("siSunPosition:input","no longitude specified in struct position"); }
	tmpVar = mxGetField(pos, 0, "altitude");
	if(tmpVar) { sData.elevation = mxGetScalar(tmpVar); } else { mexErrMsgIdAndTxt("siSunPosition:input","no altitude specified in struct position"); }
	typeField = 2;
} else {
	sData.latitude = mxGetScalar(prhs[1]);//arg2 or arg2.latitude
	sData.longitude = mxGetScalar(prhs[2]);//arg3 aka arg2.longitude
	sData.elevation = mxGetScalar(prhs[3]);//arg4 aka arg2.altitude
	typeField = 4;
}
if(nlhs > 3) {
	sData.function = SPA_ZA_RTS;
}
if(nrhs > typeField) {//We have a calculation type field
	if(mxIsChar(prhs[typeField])) {
		char * typeV = mxArrayToString(prhs[typeField]);
		if(strcasecmp(typeV,"SPA_ZA_RTS")==0) {
			sData.function = SPA_ZA_RTS;
		} else if (strcasecmp(typeV,"SPA_ZA_INC")==0) {
			sData.function = SPA_ZA_INC;
		} else if (strcasecmp(typeV,"SPA_ALL")==0) {
			sData.function = SPA_ALL;
		}
		mxFree(typeV);
	} else if(mxIsNumeric(prhs[typeField])) {
		sData.timezone = (int)mxGetScalar(prhs[typeField]);
	}
}

// Loop over time points and get data for those
m = mxGetM(prhs[0]);
n = mxGetN(prhs[0]);
// first create the output arrays
#ifdef STRUCT_OUT
plhs[0] = mxCreateStructMatrix(m, n, 4, ofns);
T_o = mxGetPr(times_o[0]);
#else
//Allocate output matrices
plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
plhs[3] = mxCreateDoubleMatrix(m, n, mxREAL);
plhs[4] = mxCreateDoubleMatrix(m, n, mxREAL);
plhs[5] = mxCreateDoubleMatrix(m, n, mxREAL);
//Get pointers so we can fill in the data
zenith = mxGetPr(plhs[0]);
azimuth = mxGetPr(plhs[1]);
es_dist = mxGetPr(plhs[2]);
sunrise = mxGetPr(plhs[3]);
sunset = mxGetPr(plhs[4]);
suntransit = mxGetPr(plhs[5]);
#endif
n*=m;

mon = t+n; day = mon+n; hr = day+n; min = hr+n; sec = min+n; //get indexes into the datevec array

for(i=0; i<n; i++) {
	#ifdef STRUCT_OUT
	mZen = mxCreateDoubleMatrix(1,1,mxREAL);
	mAz = mxCreateDoubleMatrix(1,1,mxREAL);
	mR = mxCreateDoubleMatrix(1,1,mxREAL);
	mT = mxCreateDoubleMatrix(1,1,mxREAL);
	T_ = mxGetPr(mT);
	*T_ = T_o[i];
	#endif
	sData.year = *(t++);
	sData.month = *(mon++);
	sData.day = *(day++);
	sData.hour = *(hr++);
	sData.minute = *(min++);
	sData.second = *(sec++);
	//mexWarnMsgIdAndTxt("siSP:time","Date is %04d-%02d-%02d %02d:%02d:%02d",sData.year,sData.month,sData.day,sData.hour,sData.minute,sData.second);
	spa_calculate(&sData);
	#ifdef STRUCT_OUT
	zenith = mxGetPr(mZen);
	azimuth = mxGetPr(mAz);
	es_dist = mxGetPr(mR);
	*zenith = sData.zenith;
	*azimuth = sData.azimuth;
	*es_dist = sData.r;
	mxSetFieldByNumber(plhs[0],i,0,mT);
	mxSetFieldByNumber(plhs[0],i,1,mZen);
	mxSetFieldByNumber(plhs[0],i,2,mAz);
	mxSetFieldByNumber(plhs[0],i,3,mR);
	#else
	zenith[i] = sData.zenith;
	azimuth[i] = sData.azimuth;
	es_dist[i] = sData.r;
	sunrise[i] = sData.sunrise;
	sunset[i] = sData.sunset;
	suntransit[i] = sData.suntransit;
	#endif
}

}

