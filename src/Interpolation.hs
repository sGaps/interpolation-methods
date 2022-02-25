{-|
Module      : Interpolation
Description : Common reexports
Copyright   : (c) Gabriel Peraza, 2022
License     : MIT
Stability   : experimental

This module just reexport some useful functions.
-}
module Interpolation (
    dividedDifferences,

    makeSpline,
    getSegments,
    evaluateSpline
) where

import Interpolation.CubicSpline        (makeSpline,getSegments,evaluateSpline)
import Interpolation.DividedDifferences (dividedDifferences)


