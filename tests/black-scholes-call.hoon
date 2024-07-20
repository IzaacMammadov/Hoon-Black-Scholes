/+  *test
/=  black-scholes-call  /gen/black-scholes-call
|%
++  in-between-true-prices
::  Checks whether calculated-price is in between true-price-low
::  and true-price-high
|=  [true-price-low=@rs true-price-high=@rs calculated-price=@rs]
^-  vase
!>
?&  (lth true-price-low calculated-price)
(lth calculated-price true-price-high)
==
::
++  test-01
::  Test that the price of a call option comes to approximately correct value
::
=/  true-price-low  .10.019
=/  true-price-high  .10.021
=/  calculated-price  +:(black-scholes-call [~ ~[.100 .110 .0.05 .1 .0.30] ~])
%-  expect
(in-between-true-prices true-price-low true-price-high calculated-price)
::
++  test-02
::  Test that a deeply out of the money call option comes to the
::  approximately correct value
::
=/  true-price-low  .8.229
=/  true-price-high  .8.231
=/  calculated-price  +:(black-scholes-call [~ ~[.100 .150 .0.05 .1 .0.50] ~])
%-  expect
(in-between-true-prices true-price-low true-price-high calculated-price)
::
++  test-03
::  Test that a deeply in the money call option comes to the
::  approximately correct value
::
=/  true-price-low  .34.838
=/  true-price-high  .34.840
=/  calculated-price  +:(black-scholes-call [~ ~[.100 .75 .0.05 .1 .0.50] ~])
%-  expect
(in-between-true-prices true-price-low true-price-high calculated-price)
--