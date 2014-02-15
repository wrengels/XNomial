\name{xmonte}
\alias{xmonte}
\title{Perform Multinomial Goodness-Of-Fit Test By Monte-Carlo Simulations}
\usage{
xmonte(obs, expr, ntrials = 1e+05, statName = "LLR", histobins = F,
  histobounds = c(0, 0), showCurve = T, detail = 1, safety = 1e+08)
}
\arguments{
  \item{ntrials}{the number of random trials to look at,
  such as \code{ntrials=100000}}

  \item{obs}{vector containing the observed numbers. All
  are non-negative integers summing to \code{> 0}.}

  \item{expr}{vector containing expectation. The length
  should be the same as that of \code{obs} and they should
  be non-negative summing to \code{> 0}. They need not be
  integers or sum to one.}

  \item{statName}{name of the test statistic to use as a
  measure of how deviant an observation is from the
  expectation. The choices are: \dQuote{LLR} for the
  log-likelihood ratio, \dQuote{Prob} for the probability,
  \dQuote{Chisq} for the chisquare statistic.}

  \item{histobins}{specifies histogram plot. If set to 0,
  \code{F} or \code{FALSE} no histogram is plotted. If set
  to 1 or \code{T} or \code{TRUE} a histogram with 500 bins
  will be plotted. If set to a number \code{> 1} a
  histogram with that number of bins is plotted.}

  \item{histobounds}{vector of length 2 indicating the
  bounds for the histogram, if any. If unspecified, bounds
  will be determined to include about 99.9% of the
  distribution.}

  \item{showCurve}{should an asymptotic curve be drawn over
  the histogram?}

  \item{detail}{how much detail should be reported
  concerning the P value. If 0, nothing is printed for
  cases where the function is used programmatically.
  Minimal information is printed if \code{detail} is set to
  1, and additional information if it is set to 2.}

  \item{safety}{a large number, such as one billion, to set
  a limit on how many samples will be examined. This limit
  is there to avoid long computations.}
}
\value{
\code{xmonte} returns a list with the following components:
\item{$ obs}{the observed numbers used as imput} \item{$
expr}{expected ratios} \item{$ ntrials}{the number of
random tables examined} \item{$ statType}{which test
statistic was used} \item{$ pLLR/pProb/pChi}{the P value
computed for the given test statistic} \item{$
standard.error}{the binomial standard error of the
estimated P value} \item{$ observedLLR}{the value of LLR
statistic for these data} \item{$ observedProb}{the
multinomial probability of the observed data under the null
hypothesis} \item{$ observedChi}{observed value of the chi
square statistic} \item{$ histobins}{number of bins in the
histogram (suppressed if zero)} \item{$ histobounds}{range
in histogram (suppressed if not used)} \item{$
histoData}{data for histogram (suppressed if not used)
Length is \code{histobins}} \item{$ asymptotoc.p.value}{the
P value obtained from the classical asymptotic test -- use
for comparison only}
}
\description{
Use \code{xmonte} to compute a P value to test whether a
set of counts fits a specific multinomial distribution. It
does this by examining a large number of random outcomes
and finding the probability of those cases which deviate
from the expectation by at least as much as the observed.
}
\examples{
#One of Gregor Mendel's crosses produced four types of pea seeds in the numbers:
peas <- c(315, 108, 101, 32)
#and he expected them to appear in the ratio of 9:3:3:1 according to his genetic model.
expected <- c(9, 3, 3, 1)
#Test Mendels theory using
xmonte(peas, expected)
#To see a histogram of the likelihood ratio statistic, use:
xmonte(peas, expected, histobins = TRUE)
#The red areas of the histogram represent those outcomes deviating from the expected 9:3:3:1 ratio
#at least as much as the observed numbers. (Much has been made of the tendency for Mendel's data
#to fit the expectations better than expected!)
#If you wish to use the standard chisquare statistic as a measure of goodness-of-fit instead
#of the LLR, use:
xmonte(peas, expected, statName="Chisq", histobins=TRUE)
}

