#include "modeller.h"
#include "truncated_gaussian.h"
#include <math.h>
#include <stdlib.h>

static const float rt = 0.5900991; /* RT at 297.15K, in kcal/mol */
static const float spi2 = 2.5066282731; /* sqrtf(2*3.14159265); */
static const float fact = 1.0;
static const float smallest  = 3.4e-38;
static const float largest   = 3.4e+38;

/* Decode parameters from Modeller parameter array */
static void get_param(const float *pcsr, float *w_i, float *mean, float *stdev, int modal)
{
  int iw, imean, istdv, k;
  iw = 2;
  imean = 2 + modal;
  istdv = 2 + modal + modal;
  for (k = 0; k < modal; k++) {
    iw = iw + 1;
    imean = imean + 1;
    istdv = istdv + 1;
    w_i[k] = pcsr[iw];
    mean[k] = pcsr[imean];
    stdev[k] = pcsr[istdv];
  }
}

/* Evaluate the truncated multiGaussian form */
static int myform_eval(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, gboolean deriv, float *fderv, float *val)
{
  float p_tot, ntmps, delE, delEmax, slope, xmaxscale, g_sum, g_temp,
    g_sum_max;
  float tmps[*modal + *modal + 1], savtanh[*modal + *modal + 1], gmin[*modal + 1],
    w_i[*modal + 1], mean[*modal + 1], stdev[*modal + 1], delxmax[*modal + 1],
    rv[*modal + 1], amp[*modal + 1], sumg_notk[*modal + 1], xtest[*modal + 1],
    g_toadd[*modal + 1], g_max[*modal + 1], tempF[*modal + 1];
  int k, l, m, ierr, iexpa, itest, ntest, cycle, ncycle, test_max_is_amp;

  /* Define params for PDF */
  delE = pcsr[0];
  slope = pcsr[1];
  xmaxscale = pcsr[2];
  get_param(pcsr, w_i, mean, stdev, *modal);

  p_tot = 0;
  //if assume p_tot=1 than this is not necessary
  /*for (k = 3; k < 3 + *modal; k++) { 
    p_tot = p_tot + pcsr[k];
    }*/

  ierr = 0;
  tmps[0] = 0.0;
  iexpa = 0;

  //Determine initial params for the individual Gaussians
  for (k = 0; k < *modal; k++) {
    iexpa = iexpa + 1;
    //calculate the normalized violation of the k-th basis Gaussian:
    rv[k] = mod_feature_delta(feat[0], mean[k], iftyp[0], &ierr) / stdev[k];
    if (ierr != 0) {
      return 1;
    }

    tmps[iexpa] = rv[k];
    iexpa = iexpa + 1;

    amp[k] = w_i[k] / (spi2 * stdev[k]);
    tmps[iexpa] = amp[k] * expf(-1 * (0.5 * rv[k] * rv[k]));
    //get estimate for gmin
    delEmax = delE - fact * rt * logf(amp[k]);
    gmin[k] = expf(-1 * delEmax / (fact * rt)) * w_i[k];
    delxmax[k] = stdev[k] * sqrtf(-2 * logf(gmin[k] / amp[k]));
  }

  //Determine gmin and delEmax, determine points that could be maxima to test (xtest[])
  ncycle = 20;
  for (cycle = 0; cycle < ncycle; cycle++) {
  //find xtest
  ntest = 0;
  for (k = 0; k < *modal; k++) {
    for (l = 0; l < *modal - 1; l++) {
      xtest[ntest] = (mean[k] + mean[l]) / 2.0;
      ntest = ntest + 1;
      for (m = l + 1; m < *modal; m++) {
	xtest[ntest] = (mean[k] + ((mean[l] + mean[m]) / 2.0) ) / 2.0;
	ntest = ntest + 1;
      }
    }
  }
  //find g_sum_max by testing g_sum for each xtest
  g_sum_max = 0.0;
  for (k = 0; k < *modal; k++) { g_max[k] = 0.0; }
  for (itest = 0; itest < ntest; itest++) {
    //printf("IIII %d %d %9.5f %9.5f %d\n",cycle,k,feat[0],xtest[itest],itest);
    g_sum = 0.0;
    for (k = 0; k < *modal; k++) {
      g_temp = amp[k] * expf(-0.5 * ((xtest[itest] - mean[k]) * (xtest[itest] - mean[k])) / (stdev[k] * stdev[k]) );
      g_sum = g_sum + fmaxf(g_temp, gmin[k]);
    }
    if(g_sum > g_sum_max){
      g_sum_max = g_sum;
      for (k = 0; k < *modal; k++) {
	g_temp = amp[k] * expf(-0.5 * ((xtest[itest] - mean[k]) * (xtest[itest] - mean[k])) / (stdev[k] * stdev[k]) );
	g_max[k] = fmaxf(g_temp, gmin[k]);
	tempF[k] = xtest[itest];
      }
    }
  }
  //test if max is amp[k]... **test here to make sure max_amp is used
  g_sum_max = 0;
  for (k = 0; k < *modal; k++) { 
    g_sum_max = g_sum_max + (g_max[k] - gmin[k]);
  }
  test_max_is_amp = 1;
  for (k = 0; k < *modal; k++) {
    //printf("TTTT %d %d %d %9.5f %9.5f %9.5f %9.5f %9.5f\n",cycle,k,l,feat[0],g_sum_max,g_max[k],amp[k],gmin[k]);
    if( g_sum_max > amp[k] - gmin[k] ){
      test_max_is_amp = 0; //max is not amp
    }
  }
  //get sumg_notk
  for (k = 0; k < *modal; k++) {
    sumg_notk[k] = 0.0;
    for (l = 0; l < *modal; l++) {
      if(l != k){
	if( test_max_is_amp == 1){
	  sumg_notk[k] = sumg_notk[k] + gmin[l];
	  g_toadd[k] = amp[k];
	  //printf("PPPP %d %d %d %9.5f %d\n",cycle,k,l,feat[0],test_max_is_amp);
	}
	else{
	  sumg_notk[k] = sumg_notk[k] + g_max[l];
	  g_toadd[k] = g_max[k];
	  //printf("QQQQ %d %d %d %9.5f %d\n",cycle,k,l,feat[0],test_max_is_amp);
	}
      }
    }
    //find gmin, delxmax
    //printf("OOOO %d %d %9.5f %9.5f %9.5f %9.5f %9.5f\n",itest,k,feat[0],g_toadd[k],sumg_notk[k],gmin[k],amp[k]);
    delEmax = delE - fact * rt * logf(g_toadd[k] + sumg_notk[k]); // adjusted max value of energy
    //gmin[k] = expf(-1 * delEmax / (fact * rt)) * (w_i[k] / p_tot); // corresponding cutoff for gauss
    gmin[k] = expf(-1 * delEmax / (fact * rt)) * w_i[k]; //assume p_tot=1
    delxmax[k] = stdev[k] * sqrtf(-2 * logf(gmin[k] / amp[k])); // distance between average and where gauss=gmin
  }
  }
  //k=0;
  //printf("OOOO %d %d %9.5f %9.5f %9.5f %9.5f %9.5f\n",itest,k,feat[0],g_toadd[k],sumg_notk[k],gmin[k],tempF[k]);

  //Calculate energies and forces for the truncated, multi-Gaussian
  iexpa = 0;
  for (k = 0; k < *modal; k++) {
    iexpa = iexpa + 2;

    if(feat[0] < mean[k] - delxmax[k] * xmaxscale){
      savtanh[iexpa - 1] = tanhf(slope * feat[0] + slope * (delxmax[k] - mean[k]));
      ntmps = (1 - 0.5 * (1 + savtanh[iexpa - 1])) * gmin[k] + 
	(0.5 * (1 + savtanh[iexpa - 1])) * tmps[iexpa];
    }
    else if(feat[0] > mean[k] + delxmax[k] * xmaxscale) {
      savtanh[iexpa] = tanhf(slope * feat[0] - slope * (delxmax[k] + mean[k]));
      ntmps = (1 - 0.5 * (1 + savtanh[iexpa])) * tmps[iexpa] + 
	(0.5 * (1 + savtanh[iexpa])) * gmin[k];
    }
    else{
      ntmps = tmps[iexpa];
    }

    //printf("GGGG %d %d %9.5f %9.5f\n",itest,k,feat[0],ntmps);
    tmps[0] = tmps[0] + ntmps;

  } /* end loop over basis Gaussians */

  *val = -1 * (fact * rt * logf(tmps[0]));

  if (deriv) {
    fderv[0] = 0.0;
    iexpa = -1;
    for (k = 0; k < *modal; k++) {
      iexpa = iexpa + 2;

      if(feat[0] < mean[k] - delxmax[k] * xmaxscale){
	fderv[0] = fderv[0] + 0.5 * (1 - savtanh[iexpa] * savtanh[iexpa]) * slope * (gmin[k] - tmps[iexpa + 1]) + 
	  0.5 * (1 + savtanh[iexpa]) * (tmps[iexpa + 1] * tmps[iexpa]) / stdev[k];
      }
      else if(feat[0] > mean[k] + delxmax[k] * xmaxscale) {
	fderv[0] = fderv[0] + 0.5 * (1 - savtanh[iexpa + 1] * savtanh[iexpa + 1]) * slope * (tmps[iexpa + 1] - gmin[k]) + 
	  (1 - 0.5 * (1 + savtanh[iexpa + 1])) * (tmps[iexpa + 1] * tmps[iexpa]) / stdev[k];
      }
      else{
	//p_i * rel_viol / (stdev): derivative of gauss
	fderv[0] = fderv[0] + (tmps[iexpa + 1] * tmps[iexpa]) / (stdev[k]);
      }

    }
    // const * f'(x) / f(x)
    fderv[0] = fact * rt * fderv[0] / tmps[0];
  }
  return 0;
}

/* Find index for the minimum feature mean of the truncated gaussian form */
static int myform_ivmin(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *mean,
                          int n_pcsr, int *val)
{
  float v ,vm;
  int ierr, k, ivmin;

  vm = largest;
  ivmin = 0;
  for (k = 0; k < *modal; k++) {
    v = mod_feature_delta(feat[0], mean[k], iftyp[0], &ierr);
    if(fabsf(v) < fabsf(vm)){
      vm = v;
      ivmin = k;
    }
  }
  *val = ivmin;
  return 0;
}

/* Find index for the heavy feature mean of the truncated gaussian form */
static int myform_ivheavy(void *data, const float *feat, const int *iftyp,
                            const int *modal, int n_feat, const float *w_i,
                            int n_pcsr, int *val)
{
  float w ,wm;
  int ivheav, k;

  wm = -largest;
  ivheav = 0;
  for (k = 0; k < *modal; k++) {
    w = w_i[k];
    if(fabsf(w) > fabsf(wm)){
      wm = w;
      ivheav = k;
    }
  }
  *val = ivheav;
  return 0;
}

/* Evaluate the minimum feature mean of the truncated gaussian form */
static int myform_minmean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{
  float w_i[*modal], mean[*modal], stdev[*modal];
  int ivmin;

  get_param(pcsr, w_i, mean, stdev, *modal);
  myform_ivmin(data, feat, iftyp, modal, n_feat, mean, n_pcsr, &ivmin);
  *val = mean[ivmin];
  return 0;
}

/* Evaluate the heavy feature mean of the truncated gaussian form */
static int myform_heavymean(void *data, const float *feat, const int *iftyp,
			    const int *modal, int n_feat, const float *pcsr,
			    int n_pcsr, float *val)
{
  float w_i[*modal], mean[*modal], stdev[*modal];
  int ivheav;

  get_param(pcsr, w_i, mean, stdev, *modal);
  myform_ivheavy(data, feat, iftyp, modal, n_feat, w_i, n_pcsr, &ivheav);
  *val = mean[ivheav];
  return 0;
}

/* Evaluate the minimum violation of the truncated gaussian form */
static int myform_vmin(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  float w_i[*modal], mean[*modal], stdev[*modal];
  int ierr, ivmin;
  get_param(pcsr, w_i, mean, stdev, *modal);
  myform_ivmin(data, feat, iftyp, modal, n_feat, mean, n_pcsr, &ivmin);
  *val = mod_feature_delta(feat[0], mean[ivmin], iftyp[0], &ierr);
  return ierr;
}

/* Evaluate the heavy violation of the truncated gaussian form */
static int myform_vheavy(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  float w_i[*modal], mean[*modal], stdev[*modal];
  int ierr, ivheav;
  get_param(pcsr, w_i, mean, stdev, *modal);
  myform_ivheavy(data, feat, iftyp, modal, n_feat, w_i, n_pcsr, &ivheav);
  *val = mod_feature_delta(feat[0], mean[ivheav], iftyp[0], &ierr);
  return ierr;
}

/* Evaluate the relative minimum violation of the truncated gaussian form */
static int myform_rvmin(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  float w_i[*modal], mean[*modal], stdev[*modal];
  int ierr, ivmin;
  get_param(pcsr, w_i, mean, stdev, *modal);
  myform_ivmin(data, feat, iftyp, modal, n_feat, mean, n_pcsr, &ivmin);
  *val = mod_feature_delta(feat[0], mean[ivmin], iftyp[0], &ierr) / stdev[ivmin];
  return ierr;
}

/* Evaluate the relative heavy violation of the truncated gaussian form */
static int myform_rvheavy(void *data, const float *feat, const int *iftyp,
			const int *modal, int n_feat, const float *pcsr,
			int n_pcsr, float *val)
{
  float w_i[*modal], mean[*modal], stdev[*modal];
  int ierr, ivheav;
  get_param(pcsr, w_i, mean, stdev, *modal);
  myform_ivheavy(data, feat, iftyp, modal, n_feat, w_i, n_pcsr, &ivheav);
  *val = mod_feature_delta(feat[0], mean[ivheav], iftyp[0], &ierr) / stdev[ivheav];
  return ierr;
}

/* Get the range (for splining) of the truncated gaussian form */
static int myform_range(void *data, int iftyp, int modal, const float *pcsr,
			int n_pcsr, float spline_range, float *minfeat,
                        float *maxfeat)
{
  float w_i[modal], mean[modal], stdev[modal];
  float min_mean, max_mean;
  int imin, imax, i;

  get_param(pcsr, w_i, mean, stdev, modal);

  /* Get indices of biggest and smallest means */
  imin = imax = 0;
  min_mean = max_mean = mean[0];
  for (i = 1; i < modal; ++i) {
    if (mean[i] > max_mean) {
      max_mean = mean[i];
      imax = i;
    }
    if (mean[i] < min_mean) {
      min_mean = mean[i];
      imin = i;
    }
  }

  *minfeat = mean[imin] - spline_range * stdev[imin];
  *maxfeat = mean[imax] + spline_range * stdev[imax];
  /* Assume feature is a distance, so min range cannot be negative */
  if (*minfeat < 0.) {
    *minfeat = 0.;
  }
  return 0;
}

/* Create the new truncated gaussian form, and return its identifier */
int truncated_gaussian_create(void)
{
  return mod_user_form_new2(myform_eval, NULL, myform_vmin, NULL, myform_vheavy,
                            NULL, myform_rvmin, NULL, myform_rvheavy, NULL,
                            myform_minmean, NULL, myform_heavymean, NULL,
                            myform_range, NULL);
}
