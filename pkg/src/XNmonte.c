//
//  RCmonte.c
//  RCmaker
//
//  Created by Bill Engels on 12/31/13.
//  Copyright (c) 2013 Bill Engels. All rights reserved.
//

#ifdef NOT_READY_FOR_R
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <time.h>
    #include "XNenumerated.h"
    #include "XNmonte.h"

    #define R_pow_di pow  // The following defines convert R calls to normal C
    #define Calloc(x,y) calloc(x,sizeof(y))
    #define Free free
    #define Rprintf printf
    #define lgammafn lgamma
    #define GetRNGstate() ;
    #define PutRNGstate() ;

void rmultinom(int n, double* prob, int K, int* rN) {
    int i, j;
    double x, cp;
    for (i = 0; i < K; i++) rN[i] = 0;
    for (j = 0; j < n ; j++) {
        x = (double)rand()/RAND_MAX;
        i = 0;
        cp=prob[0];
        while (x > cp) {
            cp += prob[++i];
        }
        (rN[i])++;
    }
}


#else
    #include <R.h>
    #include <Rmath.h>
#endif

double lmultiProb (int * k, double * lp, int m) {

    double x = 0;
    for (int i = 0; i < m; i++) {
        x += k[i] * lp[i] - lgammafn(1. + k[i]);
    }
    return x;
}


void montenomialTest  (int * obs,
                       double * expr,
                       int * ntrials,
                       int * nn,
                       int * statTypeR,
                       double * pLLR, // the LLR p-value
                       double * pProb, // the prob p-value
                       double * pChi, // the chi sq p-value
                       double * obsLLR, // the observed LLR
                       double * obsProb, // the observed prob
                       double * obsChiStat, // observed Chi Sq statistic
                       int * histoBinsR,
                       double * histoBounds,
                       int * histoData) {
    
    double * probs = Calloc(*nn, double);
    double * lprobs = Calloc(*nn, double);
    double gnp1; // lgamma(n+1)
    double * expected = Calloc(*nn, double);
    double statLeft, statSpan; // for histogram
    int hdex; // for histogram

    // get the total sample size, n
    unsigned n = 0;
    for (int i = 0; i < (*nn); i++) n += obs[i];
    
    // scale the exp array so that they are probabilities. (This may already have been done in R)
    double exum = 0;
    for (int i = 0; i < (*nn); i++) exum += expr[i];

    for (int i = 0; i < (*nn); i++) {
        probs[i] = expr[i]/exum;
        lprobs[i] = log(probs[i]);
        expected[i] = probs[i] * n;
    }
    
    
#ifdef NOT_READY_FOR_R
    srand((unsigned)time(NULL)); // seed machine random
    
    // compute observed values. This is normally done in R
    *obsLLR = 0;
    for (int i = 0; i < *nn; i++) {
        if (obs[i] > 0) {
            (*obsLLR) += obs[i] * log(expected[i]/obs[i]);
        }
    }
    *obsProb = exp(lmultiProb(obs, probs, *nn) + lgamma(1. + n));
    *obsChiStat = 0;
    for (int i = 0; i < *nn; i++) {
        *obsChiStat += R_pow_di(expected[i] - obs[i], 2)/expected[i];
    }
#endif
    
    
    // Adjust the observed to avoid tests for floating equality
    double adj = 1.0000000001;
    *obsProb *= adj;
    *obsLLR /= adj;
    *obsChiStat /= adj;
    
    gnp1 = lgammafn(1. + n);
    unsigned * rm = Calloc(*nn, unsigned); // Where we'll put the random multinomial
    double pr, stat;
    double lobsProb = 0;
    for (int i = 0; i < *nn; i++) {
        lobsProb += obs[i] * lprobs[i] - lgammafn(1. + obs[i]);
    }
    pr = (gnp1 + lobsProb);
    pr = exp(pr);
    lobsProb /= adj;
    double logProbPerfect = 0;
    int intexpi;
    for (int i = 0; i < *nn; i++) {
        intexpi = round(expected[i]);
        logProbPerfect += intexpi * lprobs[i] - lgammafn(1. + intexpi);
    }
    
    if (*histoBinsR) { // prepare for histogram
        for (int i = 0; i < *histoBinsR; i++) histoData[i] = 0;
        statLeft = histoBounds[0];
        statSpan = (histoBounds[1] - statLeft)/(*histoBinsR);
        if (statSpan == 0) *histoBinsR = 0; // No histogram can be made
    }
    
    *pLLR = 0; *pChi = 0; *pProb = 0;
//    GetRNGstate();

    //************************************
    //  This is the main loop to generate (*ntrials) random cases
    //************************************
    for (int kk = 0; kk < *ntrials; kk++) {
        // Get a random multinomial
        
        rmultinom(n, probs, *nn, (int*)rm);

//        // Display the random sample
//        Rprintf("\nTrial %d: ",kk);
//        for (int m = 0; m < *nn; m++) {
//            Rprintf("%5d", rm[m]);
//        }

        
        // Use switch to compute only the requested statistic to save time.
        // Actually, though, LLR and Chisquare are fast relative to getting the random multinomial. Only Prob is slow.        
        switch (*statTypeR) {
            case 1:
                stat = 0; // Use LLR as measure of distance
                for (int i = 0; i < *nn; i++) {
                    if (rm[i] > 0) {
                        stat += rm[i] * log(expected[i]/rm[i]);
                    }
                }
                if (stat <= *obsLLR) {
                    *pLLR += 1;
                }
                break;
            case 2: // Use probability of outcome as measure of "distance"
                stat = 0;
                for (int i = 0; i < *nn; i++) {
                    stat += rm[i] * lprobs[i] - lgammafn(1. + rm[i]);
                }
                if (stat <= lobsProb) {
                    *pProb += 1;
                }
                break;
            case 3:
                stat = 0; // Use chisquare as measure of distance
                for (int i = 0; i < *nn; i++) {
                    stat += R_pow_di(expected[i] - rm[i], 2)/expected[i];
                }
                if (stat >= *obsChiStat) {
                    *pChi += 1;
                }
                break;
            default:
                break;
        }
        
        if (*histoBinsR) {   //  Do this only if user requested histobram by setting *histoBinsR > 0
            if(*statTypeR == 1) stat *= (-2.); // convert to have asymptotic chisquare dist'n
            if(*statTypeR == 2) stat = -2 * (stat - logProbPerfect);
            hdex = (stat - statLeft)/statSpan;
            if ((hdex >= 0) && (hdex < *histoBinsR)) {
                (histoData[hdex])++;
            }
        }
    }

    *pLLR /= *ntrials;
    *pProb /= *ntrials;
    *pChi /= *ntrials;
    
//    PutRNGstate();
    Free(probs); Free(expected);
    Free(lprobs);
    Free(rm);
}

