import numpy as np
from abc import ABC, abstractmethod
from math import gamma

class dist(ABC):
    """ 
    A simple distribution class.
    Requires the implementation of properties `a` and `b`, representing the 
    bounds of the distribution, and methods `pdf`, `cdf` for the probability 
    and cumulative density functions and `ppf` for the inverse of `cdf`.
    """

    # @property
    # @abstractmethod
    # def a(self): pass
    
    # @property
    # @abstractmethod
    # def b(self): pass

    @abstractmethod
    def pdf(self, x): 
        """ Probability density function evaluated at x. """
        pass

    @abstractmethod
    def cdf(self, x): 
        """ Cumulative density function evaluated at x. """
        pass

    @abstractmethod
    def ppf(self, x): 
        """ Percent point function (inverse of `cdf`) evaluated at q. """
        pass

    def sf(self, x):
        """ Survival function (1 - `cdf`) evaluated at x. """
        return 1 - self.cdf(x)

    def rvs(self, size=1):
        """
        Random samples from the distribution.

        Parameters
        ----------
        size : int or tuple of ints, optional
            Number of random samples to draw (default is 1).

        Returns
        -------
        rvs : ndarray or scalar
        """
        u = np.random.uniform(size=size)
        s = self.ppf(u)
        return np.asscalar(s) if size==1 else s

    def interval(self, alpha):
        """
        Confidence interval with equal areas around the median.

        Parameters
        ----------
        alpha : array_like of float
            Probability that an rv will be drawn from the returned range.
            Each value should be in the range [0, 1].

        Returns
        -------
        a, b : ndarray of float
            end-points of range that contains 100 * alpha % of the
            distribution support.
        """

        alpha = np.asarray(alpha)
        if np.any((alpha > 1) | (alpha < 0)):
            raise ValueError("alpha must be between 0 and 1 inclusive")
        q1 = (1.0-alpha)/2
        q2 = (1.0+alpha)/2
        a = self.ppf(q1)
        b = self.ppf(q2)
        return a, b


class kumaraswamy(dist):
    """ 
    A Kumaraswamy distribution, similar in many ways to the Beta distribution.

    This distribution is defined over the (0, 1) interval, and has parameters
    a and b. The probability density function is
  
    ```
    pdf(x; a, b) = a * b * x**(a - 1) * (1 - x**a)**(b -  1)
    ```

    Parameters
    ----------
    
    Methods
    -------
    pdf    : probability density function
    logpdf : logarithm of the probability density function
    cdf    : cumulative density function
    ppf    : percent point function (inverse of cdf)

    Properties
    ----------
    mean, mode, var, std, skewness, kurtosis

    """

    def __init__(self, a, b):
        assert a>0 and b>0, \
            'parameters `a` and `b` must both be positive'
        self.a, self.b = a, b

    def _support_mask(self, x):
        return (0 <= x) & (x <= 1)
    def _separate_support_mask(self, x):
        return 0 <= x, x <= 1

    def pdf(self, x):
        x = np.asarray(x)
        m = self._support_mask(x)
        a, b = self.a, self.b
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(m, a * b * x**(a-1) * (1-x**a)**(b-1), 0.)
    
    def logpdf(self, x):
        pass
        x = np.asarray(x)
        m = self._support_mask(x)
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.log(self.a) + np.log(self.b) + (self.a-1)*np.log(x) + (self.b-1)*np.log(1-x**self.a)

    def cdf(self, x):
        x = np.asarray(x)
        m1, m2 = self._separate_support_mask(x)
        with np.errstate(invalid='ignore'):
            return np.where(m2, np.where(m1, 1 - (1 - x**self.a)**self.b, 0.), 1.)
    
    def ppf(self, p):
        p = np.asarray(p)
        m1, m2 = 0 < p, p < 1
        return np.where(m2, np.where(m1, (1 - (1-p)**(1/self.b))**(1/self.a), 0.), 1.)

    @property
    def mean(self):
        return (self.b * gamma(1+1/self.a) * gamma(self.b)) / gamma(1 + 1/self.a + self.b)
    
    @property
    def mode(self):
        """
        Note that
        - if a>1 and b>1 the distribution is unimodal
        - if a<1 and b<1 the distribution is uniantimodal
        - if a>1 and b<=1 the distribution is increasing
        - if a<=1 and b>1 the distribution is decreasing
        - if a=b=1 the distribution is constant

        This property is only defined in the first two cases, returning the 
        mode and the antimode, respectively.
        """ 
        if (self.a>1 and self.b>1) or (self.a<1 and self.b<1):
            return ((self.a - 1) / (self.a * self.b - 1))**(1 / self.a)
        else:
            return np.nan
    
    def moment(self, n):
        """ Raw moment of order n """
        return (self.b * gamma(1+n/self.a) * gamma(self.b)) / gamma(1 + n/self.a + self.b)

    @property
    def var(self):
        return self.moment(2) - self.moment(1)**2
    
    @property
    def std(self):
        return np.sqrt(self.var)
    
    @property
    def skewness(self):
        raise(NotImplementedError)
    
    @property
    def kurtosis(self):
        raise(NotImplementedError)


