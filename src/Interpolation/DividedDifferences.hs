{-|
Module      : Interpolation.DividedDifferences
Description : This module implements the divided differece method.
Copyright   : (c) Gabriel Peraza, 2022
License     : MIT
Stability   : experimental

This module implements the divided differece method. It also exports
a handy function to obtain the coefficients directly, which is
what we're interested often.

-}
module Interpolation.DividedDifferences (
    coefficients,
    dividedDifferences
) where

import Data.List (scanl', tails)

-- | Owl-Operator, combines an unary operator with a binary operator
--
-- prop> f (g x y) = (f .: g) x y
(.:) = ((.) . (.))

-- | Returns the coefficients of the divided Differences Method
--
-- See also: `dividedDifferences`
coefficients :: RealFrac a => [a] -> [a] -> [a]
coefficients = fmap head .: dividedDifferences

-- | Returns the transposed matrix of coefficents belonging to the divided differences method
--  
-- @
--
-- c00      c10 ..          c{n-1}0
-- c01      c11 ..          c{n-1}1
-- c02      c12 ..        .´
-- .        .           .´
-- .        .         .´
-- .        .       .´
-- c0{n-1}  c1{n-1}
-- cn0
--  ^
--  |
-- coefficients
--
-- @
--
-- The result is correct iff
--
-- prop> length xs == length ys
dividedDifferences :: RealFrac a => [a] -> [a] -> [[a]]
dividedDifferences ys xs = scanl' combine ys . distances $ xs
    where combine  cs = zipWith (/) (differences cs cs)
          differences = zipWith (-) . tail
          
-- | Returns a list with the following form:
-- @
-- [ [x1-x0, x1-x2, ...] , [x0-x2, x1-x3, ...] ... ]
-- @
distances :: Num a => [a] -> [[a]]
distances xs = init . tail
                    . fmap (flip (zipWith (-)) xs)
                    . tails
                    $ xs

