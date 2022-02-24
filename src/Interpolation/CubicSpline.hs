module Interpolation.CubicSpline (
    Spline (..),
    InclusiveRange (..),
    SplineCoefficients (..),
    SplineSegment,

    makeSpline,
    showSegment,
    evaluateSpline,
    evaluateSegment,

    naturalCubicSpline,
) where

import Data.List (scanl',zip3,zipWith3,zip4,find)
import Data.Function (on)
import Data.Tuple.Extra (uncurry3)

-- | Holds the coefficients (Ai, Bi, Ci, yi)
newtype SplineCoefficients a = Coefficients (a,a,a,a)
    deriving( Eq , Ord , Show )

-- Represents an inclusive Range [low, high]
newtype InclusiveRange a = IRange (a,a)
    deriving( Eq , Ord , Show )

-- | Holds The Range, The coefficient and the root for the segment's polynomial
data SplineSegment a = Segment (InclusiveRange a) (SplineCoefficients a) a
                       deriving( Eq , Ord , Show )

-- | Represents an Spline
newtype Spline a = Spline { getSegments :: [SplineSegment a] }
                   deriving( Eq , Ord , Show )

-- |
showSegment :: Show a => SplineSegment a -> String -> String
showSegment (Segment (IRange (lo,hi)) (Coefficients (a,b,c,y)) t) =
    showString "<Range: [" . shows lo . showString ", " . shows hi . showChar ']'
                        . showString ", root (t_i): " . shows t
                        . showString ", a_i: " . shows a
                        . showString ", b_i: " . shows b
                        . showString ", c_i: " . shows c
                        . showString ", y_i: " . shows y
                        . showChar '>'


-- | Bottom Value for evaluateSpline
absurdValue :: Fractional a => a
absurdValue = 1/0

-- not the best way to evaluate a spline but it works.
-- (The best representation could use a Finger Tree to order the given ranges)

-- | Evaluates x on the given Spline
--
-- See Also `makeSpline` and `naturalCubicSpline`
evaluateSpline :: (Ord a, Fractional a) => a -> Spline a -> a
evaluateSpline x = safeHead . fmap (evaluateSegment x)
                            . take 1
                            . filter (x `insideOf`)
                            . getSegments
    where insideOf x (Segment (IRange (low,high)) _ _) =
                low <= x && x <= high
          safeHead xs = if null xs then absurdValue else head xs

evaluateSegment x (Segment _ (Coefficients (ai,bi,ci,yi)) ti) =
    yi + (x - ti)*(ci + (x - ti)*(bi + (x - ti)*ai))

-- | Builds an spline with a list of ys and knots (ts)
makeSpline ys ts = Spline . zipWith toSegment (previousCurrent ts)
                          . zipWith toCoefficients ys
                          . uncurry3 zip3 
                          . naturalCubicSpline ys
                          $ ts
    where toCoefficients y (a,b,c)  = Coefficients (y,a,b,c)
          toSegment (x,xnext) coeff = Segment (IRange (x,xnext)) coeff x

---------------------
-- the actual method:

-- | Applies xi - x_{i-1} for each element of the given list
selfDifference xs = (zipWith (-) . drop 1) xs xs


-- | Returns a list of (xi, x_{i+1})
previousCurrent = zip <*> tail      -- uses s-combinator on the function zip

-- | returns hi and bi terms of the method.
prepareTridiagonal :: Fractional a => [a] -> [a] -> ([a], [a])
prepareTridiagonal ys ts =
    let hs = selfDifference ts
        bs = (zipWith merge hs . selfDifference) ys
    in (hs, bs)
    where merge h ydiff = 6 * ydiff / h

-- | returns the ui and vi terms of the method
--
-- See also `prepareTridiagonal`
tridiagonal :: Fractional a => [a] -> [a] -> ([a], [a])
tridiagonal hs bs =
    let (h0 : h1 : _) = hs
        (b0 : b1 : _) = bs
        u1            = 2 * (h0 + h1)
        v1            = b1 - b0

    in unzip . scanl' fowardReplace (u1, v1)
             . (zip `on` (tail . previousCurrent)) hs
             $ bs

    where fowardReplace (upre,vpre) ((hpre, h), (bpre, b)) =
            let u = 2 * (h + hpre) - (hpre^2)/upre
                v = b - bpre - (hpre * vpre/upre)
            in (u, v)

-- | returns the list of zi values of the cubicSpline method.
--
-- See also `tridiagonal`
backwardSustitution :: Fractional a => a -> a -> [a] -> [a] -> [a] -> [a]
backwardSustitution z0 zn hs vs us = 
    let hsBody = drop 1 hs
    in z0 : scanr sustitution zn (zip3 hsBody vs us)
    where sustitution (h,v,u) znext = v - (h*znext)/u

-- | Returns a list of coefficients. It requires to be called after
-- other methods.
--
-- See also `backwardSustitution`
coefficients :: Fractional a => [a] -> [a] -> [a] -> ([a], [a], [a])
coefficients ys hs zs =
    let as = zipWith sustitutionA hs (selfDifference zs)
        bs = fmap (/2) zs
        cs = zipWith3 sustitutionC hs (previousCurrent zs) (previousCurrent ys)
    in (as,bs,cs)
    where sustitutionA h zdiff = zdiff/(6*h)
          sustitutionC h (z,znext) (y,ynext) =
            (-h/6)*znext - (h/3)*z + (1/h)*(ynext - y)

-- | Returns the coefficients of a natural cubic spline
naturalCubicSpline :: Fractional a => [a] -> [a] -> ([a], [a], [a])
naturalCubicSpline ys ts = coefficients ys hs zs
    where (hs, bs) = prepareTridiagonal ys ts
          (us, vs) = tridiagonal hs bs
          (z0, zn) = (0,0) -- property of natural cubic splines.
          zs       = backwardSustitution z0 zn hs vs us

