:: Contains useful arms using Black-Scholes Option Pricing formulae.
|%
++  calculate-call-option-price
|=  [spot=@rs strike=@rs rate=@rs time-to-expiry=@rs implied-volatility=@rs]
::  spot: @rs - The current price in USD of the underlying asset in the
::              spot market
::  strike: @rs - The strike/exercise price of the Option contract in USD
::  rate: @rs - The continuously compounded risk-free rate per year at
::              which USD can be borrowed or lent
::  time-to-expiry: @rs - The number of years until expiry/maturity of
::                        the option contract
::  implied-volatility: @rs - The volatility of the undering asset.
::
=/  d1
%+  div:rs
%+  add:rs  (log (div:rs spot strike))
%+  mul:rs  time-to-expiry
%+  sub:rs  rate
(div:rs (mul:rs implied-volatility implied-volatility) .2)
(mul:rs implied-volatility (sqt:rs time-to-expiry))
=/  d2
(add:rs d1 (mul:rs implied-volatility (sqt:rs time-to-expiry)))
%+  sub:rs
(mul:rs (normal-cdf d2) spot)
:(mul:rs (normal-cdf d1) strike (exponential :(mul:rs .-1 rate time-to-expiry)))
::
++  normal-cdf
::  The cumulative distribution function of the Normal distribution with
::  mean 0 and variance 1. It is calculated by using the complementary
::  error function.
::
|=  x=@rs
^-  @rs
%+  mul:rs  .0.5
%-  complementary-error-function
%+  sub:rs  .0
%+  div:rs  x
%-  sqt:rs  .2
::
++  complementary-error-function
::  The complementary error function erfc(x) is defined as:
::  2/sqrt(pi) * integral x->inf e^(-t^2) dt.
::  The implementation below is sourced from ``Numerical Recipes. The Art
::  of Scientific Computing - 3rd Edition``, page 265. The algorithm has
::  a uniform error of less than 1.2e-7.
::
|=  x=@rs
^-  @rs
?.  (sig:rs x)
  (sub:rs .2 (complementary-error-function (mul:rs .-1 x)))
=/  t  (div:rs .2 (add:rs .2 x))
%+  mul:rs  t
%-  exponential
%+  add:rs
;:  mul:rs  .-1  x  x  ==
%^  fma:rs  t
%^  fma:rs  t
%^  fma:rs  t
%^  fma:rs  t
%^  fma:rs  t
%^  fma:rs  t
%^  fma:rs  t
%^  fma:rs  t
%^  fma:rs  t
.0.17087277
.-0.82215223
.1.48851587
.-1.13520398
.0.27886807
.-0.18628806
.0.09678418
.0.37409196
.1.00002368
.-1.26551223
::
++  log
::  The natural log function log(x) [sometimes written ln(x)] with base e.
::  Implemented by pulling out facors of 2, and then applying the
::  Taylor expanson for log(1+x) valid for x in (-1, 1). The Taylor expansion
::  continues until the
::
|=  x=@rs
^-  @rs
?:  (gth:rs x (div:rs .4 .3))
  (add:rs (log (div:rs x .2)) log-2)
?:  (lte:rs x (div:rs .2 .3))
  (sub:rs (log (mul:rs x .2)) log-2)
=/  z  (sub:rs x .1)
=/  n  2
=/  total  .0
=/  increment  z
|-
?:  &((lth:rs increment .1e-10) (gth:rs increment .-1e-10))
  (add:rs total increment)
%=
$
n  +(n)
total  (add:rs total increment)
increment  (mul:rs (minus-one-to-power +(n)) (div:rs (power z n) (sun:rs n)))
==
::
++  log-2  .0.69314718056
::  log(2) -- stored for general calculations of log(x)
::
++  exponential
::  The exponential function exp(x) defined for @rs. It adds terms in the
::  Taylor expansion until the term added is less than 1e-10.
::
|=  x=@rs
^-  @rs
=/  n  1
=/  total  .0
=/  increment  .1
|-
?:  &((lth:rs increment .1e-10) (gth:rs increment .-1e-10))
  (add:rs total increment)
%=
$
n  +(n)
total  (add:rs total increment)
increment  (div:rs (power x n) (sun:rs (factorial n)))
==
::
++  factorial
::  The factorial function x! defined for @ud.
::
|=  x=@ud
^-  @ud
?~  x
  1
(mul x (factorial (sub x 1)))
::
++  power
::  A @rs x raised to an integer power n. Uses O(log(n)) algorithm for
::  exponentiation.
::
|=  [x=@rs n=@ud]
^-  @rs
?~  n
  .1
?~  (mod n 2)
  (power (mul:rs x x) (div n 2))
(mul:rs x (power (mul:rs x x) (div (sub n 1) 2)))
::
++  minus-one-to-power
::  Return (-1)^x as a @rs without rounding
::
|=  x=@ud
^-  @rs
?~  (mod x 2)
  .1
.-1
--