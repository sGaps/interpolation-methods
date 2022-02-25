{-|
Module      : Interpolation.CubicSpline
Description : This module implements the Cubic Splines method.
Copyright   : (c) Gabriel Peraza, 2022
License     : MIT
Stability   : experimental

This module implements the cubic spline method. It also exports
some functions that makes easier to evaluate and print the spline.
-}
module Interpolation.CubicSpline (
    Spline(..),
    Segment(..),

    naturalSpline,
    evaluateSpline,
    evaluateSegment,
    makeSpline,

    showSegment,
    strSegment
) where

import Control.Monad.ST
import Control.Monad (forM_)
import qualified Data.Vector         as V
import qualified Data.Vector.Mutable as MV


-- Returns a vector of hs' and bs'
prepareTridiagonal :: Fractional a => Int -> V.Vector a -> V.Vector a -> (V.Vector a, V.Vector a)
prepareTridiagonal size ts ys = V.unzip $ runST $ do
    let limit = size-2
    hsbs <- MV.new size

    forM_ [0 .. limit] (\i -> do
        tnext <- V.indexM ts (i+1)
        t     <- V.indexM ts (i)

        ynext <- V.indexM ys (i+1)
        y     <- V.indexM ys (i)
        let h = tnext - t
            b = 6 * (ynext - y) / h
        MV.write hsbs (i) (h,b)
        )

    V.freeze hsbs

tridiagonal :: Fractional a => Int -> V.Vector a -> V.Vector a -> (V.Vector a, V.Vector a)
tridiagonal size hs bs = V.unzip $ runST $ do
    let limit = size-2
        
        u1 = 2*(hs V.! 0 + hs V.! 1)
        v1 = (bs V.! 1) - (bs V.! 0)

    usvs <- MV.new size

    -- sets a default value for the algorithm's 'uninitialised' element.
    MV.write usvs (0) (0,0)

    MV.write usvs (1) (u1,v1)

    forM_ [2 .. limit] (\i -> do
        hpre  <- V.indexM hs (i-1)
        h     <- V.indexM hs (i)

        bpre  <- V.indexM bs (i-1)
        b     <- V.indexM bs (i)

        (upre,vpre) <- MV.read usvs (i-1)

        let u = 2 * (h + hpre) - (hpre^2)/upre
            v = b - bpre - (hpre*vpre/upre)

        MV.write usvs (i) (u,v)
        )
    V.freeze usvs

backwardSustitution :: Fractional a => Int -> a -> a -> V.Vector a -> V.Vector a -> V.Vector a -> V.Vector a
backwardSustitution size z0 zn hs us vs = runST $ do
    let limit = size-2

    zs <- MV.new size
    MV.write zs (0) 0
    MV.write zs (size-1) 0

    forM_ [limit, limit-1 .. 1] (\i -> do
        h <- V.indexM hs (i)
        u <- V.indexM us (i)
        v <- V.indexM vs (i)

        znext <- MV.read zs (i+1)

        let z = (v - h*znext)/u
        MV.write zs (i) z
        )
    V.freeze zs

coefficientA :: Fractional a => Int -> V.Vector a -> V.Vector a -> V.Vector a
coefficientA size hs zs = runST $ do
    let limit = size-1
    as <- MV.new size

    forM_ [0 .. limit] (\i -> do
        h     <- V.indexM hs (i)
        z     <- V.indexM zs (i)
        znext <- V.indexM zs (i+1)
        let a = (1/6*h)*(znext - z)
        MV.write as (i) a
        )
    V.freeze as

coefficientB :: Fractional a => V.Vector a -> V.Vector a
coefficientB zs = fmap (/2) zs

coefficientC :: Fractional a => Int -> V.Vector a -> V.Vector a -> V.Vector a -> V.Vector a
coefficientC size ys hs zs = runST $ do
    let limit = size-1
    cs <- MV.new size

    forM_ [0 .. limit] (\i -> do
        h     <- V.indexM hs (i)

        z     <- V.indexM zs (i)
        znext <- V.indexM zs (i+1)

        y     <- V.indexM ys (i)
        ynext <- V.indexM ys (i+1)
        let c = (-h/6)*znext - (h/3)*z + (1/h)*(ynext - y)
        MV.write cs (i) c
        )
    V.freeze cs


naturalSpline :: Fractional a => [a] -> [a] -> (V.Vector a, V.Vector a, V.Vector a)
naturalSpline listTs listYs = (coefas,coefbs,coefcs)
    where n       = V.length ts
          (ts,ys) = (V.fromList listTs, V.fromList listYs)
          (hs,bs) = prepareTridiagonal n ts ys
          (us,vs) = tridiagonal n hs bs
          (z0,zn) = (0,0)

          zs      = backwardSustitution n z0 zn hs us vs
          coefas  = coefficientA (n-1) hs zs
          coefbs  = coefficientB zs
          coefcs  = coefficientC (n-1) ys hs zs

-- | Coefficients root inclusive_range and (a,b,c,y) for each S(x):
-- where S(x) = y + (x - t)(c + (x - t)(b + (x - t)a))
-- and (low,high) = inclusive_range && low <= x && x <= high
data Segment p = Segment p (p,p) (p,p,p,p)
                    deriving (Show,Eq,Ord)

-- | Type that holds the Spline's segments.
newtype Spline p = Spline { getSegments :: V.Vector (Segment p) }
                    deriving (Show,Eq)

makeSpline :: Fractional a => [a] -> [a] -> Spline a
makeSpline ts ys = Spline . V.zipWith3 Segment ts' (prevCurr ts') . V.zip4 as bs cs $ ys'
    where ys'        = V.fromList ys
          ts'        = V.fromList ts
          (as,bs,cs) = naturalSpline ts ys
          prevCurr   = V.zip <*> V.drop 1

absurdValue :: Fractional a => a
absurdValue = 1/0

evaluateSegment x (Segment t _ (a,b,c,y) ) =
    y + (x - t)*(c + (x - t)*(b + (x - t)*a))

evaluateSpline :: (Fractional a , Ord a) => a -> Spline a -> a
evaluateSpline x spline = maybe absurdValue (evaluateSegment x)
                                . V.find (inRange x)
                                . getSegments $ spline
    where inRange x (Segment _ (low,high) _) = low <= x && x <= high

showSegment :: Show a => Segment a -> String -> String
showSegment (Segment t (lo,hi) (a,b,c,y)) =
    showString "<Range: [" . shows lo . showString ", " . shows hi . showChar ']'
                        . showString ", root (t_i): " . shows t
                        . showString ", a_i: " . shows a
                        . showString ", b_i: " . shows b
                        . showString ", c_i: " . shows c
                        . showString ", y_i: " . shows y
                        . showChar '>'

strSegment :: Show a => Segment a -> String
strSegment s = showSegment s ""
