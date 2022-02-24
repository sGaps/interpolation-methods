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



-- | Owl-Operator, combines an unary operator with a binary operator
--
-- prop> f (g x y) = (f .: g) x y
(.:) = ((.) . (.))

-- vectorial sum
--2- (.+.), (.-.), (.*.) :: Num a => [a] -> [a] -> [a]
--2- (.+.) = zipWith (+)
--2- (.-.) = zipWith (-)
--2- (.*.) = zipWith (*)
--2- 
--2- (./.) :: Fractional a => [a] -> [a] -> [a]
--2- (./.) = zipWith (/)
--2- 
--2- selfDifference xs = ((.-.) . tail) xs xs

-- v 1
-- -- hs = selfDifference xs
-- -- bs = ((./.) . fmap (6*) . selfDifference) ys hs
-- -- u1 = 2 * (hs !! 0 + hs !! 1)
-- -- v1 = bs !! 1 - bs !! 0
-- -- (vs,us) = (unzip . scanl' tridiagonalUs u1) ((zip . drop 1) hs hs) ((zip . drop 1) bs bs)
-- --     where tridiagonalUs (vi_1,ui_1) ((hi,hi_1), (bi,bi_1)) =
-- --                     (2 * (hi + hi_1) - (hi_1 ^ 2)/ui_1,
-- --                      bi - bi_1 - (hi_1 * vi_1 / ui_1) )
-- -- 
-- -- -- natural spline: both extermes are zero
-- -- zn = 0
-- -- zs = 0 : scanr backward zn (zip3 hs vs us)
-- --     where backward (hi, vi, ui) zi_plus1 = (vi - hi* zi_plus1)/ui
-- -- 
-- -- coefAs = zipWith (*) (fmap (\h -> 1/(6*h)) hs) (selfDifference zs)
-- -- coefBs = fmap (/2) zs
-- -- coefCs = zs_foward .+. zs_normal .+. ys_diffs
-- --     where zs_foward = fmap (\zh -> -zh/6) ((zipWith (*) . drop 1) zs hs)
-- --           zs_normal = fmap (\zh -> -zh/3) (zipWith (*) zs hs)
-- --           ys_diffs  = zipWith (/) ((zipWith (*) . drop 1) ys ys) hs

-- v3

-- | (Ai, Bi, Ci, Di) or (Ai, Bi, Ci, yi)
--type SplineCoefficient a = (a, a, a, a)
--
---- | [low, high]
--type InclusiveRange a = (a, a)

-- | (Ai, Bi, Ci)
newtype SplineCoefficients a = Coefficients (a,a,a,a)
    deriving( Eq , Ord , Show )

-- [low, high]
newtype InclusiveRange     a = IRange (a,a)
    deriving( Eq , Ord , Show )

-- these are for the driver program
--data SplineSegment a = Segment (InclusiveRange a) (SplineCoefficient a) a
data SplineSegment a = Segment (InclusiveRange a) (SplineCoefficients a) a
                       deriving( Eq , Ord , Show )

--showSegment (Segment rng (a,b,c) y t) =
showSegment (Segment (IRange (lo,hi)) (Coefficients (a,b,c,y)) t) =
    showString "<Range: [" . shows lo . showString ", " . shows hi . showChar ']'
                        . showString ", root (t_i): " . shows t
                        . showString ", a_i: " . shows a
                        . showString ", b_i: " . shows b
                        . showString ", c_i: " . shows c
                        . showString ", y_i: " . shows y
                        . showChar '>'

newtype Spline a = Spline { getSegments :: [SplineSegment a] }
                   deriving( Eq , Ord , Show )
-- [Segment (low,high) (a,b,c,y) t, ...]

absurdValue :: Fractional a => a
absurdValue = 1/0

-- not the best way to evaluate a spline but it works.
-- (The best representation uses a Finger Tree to order the given ranges)
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

-- takes ys ts
makeSpline ys ts = Spline . zipWith toSegment (previousCurrent ts)
                          . zipWith toCoefficients ys
                          . uncurry3 zip3 
                          . naturalCubicSpline ys
                          $ ts
    where toCoefficients y (a,b,c)  = Coefficients (y,a,b,c)
          toSegment (x,xnext) coeff = Segment (IRange (x,xnext)) coeff x



--makeSpline ys ts =
--    let (as,bs,cs) = naturalCubicSpline ys ts
--        mixed      = zip3 as bs cs
--    in zipWith toCoefficients ys $ mixed
--    where toCoefficients y (a,b,c)  = Coefficients (y,a,b,c)
--          toSegment (x,xnext) coeff = Segment (IRange (x,xnext)) coeff x

    

---------------------
-- the actual method:

selfDifference xs = (zipWith (-) . drop 1) xs xs

-- uses s-combinator over functions (->)
previousCurrent   = zip <*> tail 

-- it seems that work fine...
-- | returns hi and bi terms of the method.
prepareTridiagonal ys ts =
    let hs = selfDifference ts
        bs = (zipWith merge hs . selfDifference) ys
        --bs = zipWith merge hs (previousCurrent ys)
    in (hs, bs)
    where merge h ydiff = 6 * ydiff / h
          --merge h (ypre,y) = 6 * (y - ypre) / h
          previousCurrent  = zip <*> tail -- uses s-combinator over functions (->)

-- expected b0 = (24)/ (1)
-- expected b1 = 6*(3)/ (1) = 18
-- expected b2 = 6*(2-9)/ (1) = 6*(-7) = -42

-- | returns the ui and vi terms of the method
tridiagonal hs bs =
    let (h0 : h1 : _) = hs
        (b0 : b1 : _) = bs
        u1            = 2 * (h0 + h1)
        v1            = b1 - b0

    -- TODO: it starts from h1 and b1, so it must be zip (drop 1 (combined)
    in unzip . scanl' fowardReplace (u1, v1)
             . (zip `on` (tail . previousCurrent)) hs -- TODO: Do i have to use tail here?
             $ bs

    where fowardReplace (upre,vpre) ((hpre, h), (bpre, b)) =
            let u = 2 * (h + hpre) - (hpre**2)/upre
                v = b - bpre - (hpre * vpre/upre)
            in (u, v)
          previousCurrent = zip <*> tail -- uses s-combinator over functions (->)

-- seems good until tridiagonal
step ys ts =
    let (hs , bs) = prepareTridiagonal ys ts
        (us , vs) = tridiagonal hs bs
        zs = backwardSustitution 0 0 hs vs us
    --in (zs , (us,vs) , (hs, bs))
    in (zs , (us,vs,hs) )

-- We have that:
--      length u == length v
-- but hs has two extra element, so
--      length u == length h - 2
--
-- And if that weren't enougth, this algorithm iterates from the value z1 until z_{n-1},
-- so it ignores h0 (the hs' head):
--
-- Visually, he have that:
--      us = [u1, u2, ... u_{n-1}]
--      vs = [v1, v2, ... v_{n-1}]
--
-- Which means that we have to crop hs' head:
--      hs = h0 : [h1, h2, ... , h_{n-1}
--
-- to obtain the desired result. woah!

-- | returns the list of zi values of the cubicSpline method.
backwardSustitution z0 zn hs vs us = 
    --let hsBody = (init . tail) hs
    let hsBody = drop 1 hs
    in z0 : scanr sustitution zn (zip3 hsBody vs us)
    where sustitution (h,v,u) znext = v - (h*znext)/u

coefficients :: RealFrac a => [a] -> [a] -> [a] -> ([a], [a], [a])
coefficients ys hs zs =
    --let as = zipWith sustitutionA hs (previousCurrent zs)
    let as = zipWith sustitutionA hs (selfDifference zs)
        bs = fmap (/2) zs
        cs = zipWith3 sustitutionC hs (previousCurrent zs) (previousCurrent ys)
    in (as,bs,cs)
    where previousCurrent      = zip <*> tail -- uses s-combinator over functions (->)
          sustitutionA h zdiff = zdiff/(6*h)
          --sustitutionA h (z,znext) = (znext-z)/(6*h)
          sustitutionC h (z,znext) (y,ynext) =
            (-h/6)*znext - (h/3)*z + (1/h)*(ynext - y)

-- 8 points, implies 7 equations
-- but 6 are returned... Something is odd
naturalCubicSpline ys ts = coefficients ys hs zs
    where (hs, bs) = prepareTridiagonal ys ts
          (us, vs) = tridiagonal hs bs
          (z0, zn) = (0,0) -- property of natural cubic splines.
          zs       = backwardSustitution z0 zn hs vs us


-- | It takes a list of f(ti) and a list of knots ts and returns a list of coefficients
cubicSpline'' ys ts = (coefAs,coefBs,coefCs)
    where -- tridiagonal solver function
          tridiagonal (vi_1,ui_1) ((hi, hi_1), (bi, bi_1)) =
                              let ui = 2 * (hi + hi_1) - (hi_1 ^ 2)/ui_1
                                  vi = bi - bi_1 - (hi_1 * vi_1 / ui_1)
                              in (vi, ui)

          -- backward sustitution function
          backwardSust (hi, vi, ui) zi_plus1 = (vi - hi* zi_plus1)/ui

          hs      = selfDifference ts
          bs      = (zipWith (/) . fmap (6*) . selfDifference) ys hs
          u1      = 2 * (hs !! 0 + hs !! 1)
          v1      = bs !! 1 - bs !! 0
          -- BUG: it uses (u1, v1) in `scanl'` instead of (v1,u1)
          (vs,us) = unzip . scanl' tridiagonal (v1,u1) . zip (tail hs `zip` hs) $ (tail bs `zip` bs)

          zn = 0
          z0 = 0
          zs = z0 : scanr backwardSust zn (zip3 hs vs us)

          coefAs = fmap (\(zdiff,h) -> zdiff/(6*h) ) . zip (selfDifference zs) $ hs
          coefBs = fmap (/2) zs
          
          zs_foward = fmap (\zh -> -zh/6) ((zipWith (*) . drop 1) zs hs)
          zs_normal = fmap (\zh -> -zh/3) (zipWith (*) zs hs)
          ys_diffs  = zipWith (/) ((zipWith (*) . drop 1) ys ys) hs

          -- check the zs_foward
          coefCs = zipWith (+) zs_foward . zipWith (+) zs_normal $ ys_diffs

-- v2
--2- hs      = selfDifference xs
--2- bs      = (zipWith (/) . fmap (6*) . selfDifference) ys hs
--2- u1      = 2 * (hs !! 0 + hs !! 1)
--2- v1      = bs !! 1 - bs !! 0
--2- (vs,us) = (unzip . scanl' tridiagonal u1) (tail hs `zip` hs) (tail bs `zip` bs)
--2-     where tridiagonal (vi_1,ui_1) ((hi, hi_1), (bi, bi_1)) =
--2-                     let ui = 2 * (hi + hi_1) - (hi_1 ^ 2)/ui_1
--2-                         vi = bi - bi_1 - (hi_1 * vi_1 / ui_1)
--2-                     in (vi, ui)
--2- zn = 0
--2- z0 = 0
--2- zs = z0 : scanr backwardSust zn (zip3 hs vs us)
--2-     where backwardSust (hi, vi, ui) zi_plus1 = (vi - hi* zi_plus1)/ui
--2- 
--2- coefAs = fmap (\(zdiff,h) -> zdiff/(6*h) ) (zip . selfDifference) zs hs
--2- coefBs = fmap (/2) zs
--2- -- check the zs_foward
--2- coefCs = zs_foward .+. zs_normal .+. ys_diffs
--2-     where zs_foward = fmap (\zh -> -zh/6) ((zipWith (*) . drop 1) zs hs)
--2-           zs_normal = fmap (\zh -> -zh/3) (zipWith (*) zs hs)
--2-           ys_diffs  = zipWith (/) ((zipWith (*) . drop 1) ys ys) hs


--selfDifference xs = ((.-.) . tail) xs xs
--
--hs = (zipWith (-) . drop 1) ts ts
--us = (fmap (2*) . zipWith (+) . drop 1) hs hs
--bs = flip (zipWith (/)) hs . fmap (6*) . (zipWith (-) . drop 1) ys ys
--vs = (zipWith (-) . drop 1) bs bs
--
--zs = zipWith (-) . 

-----
--differences       = zipWith (-) . tail
--selfDifference xs = differences xs xs

--setup ys ts = (hs,bs)
--    where hs = selfDifference ts
--          bs = zipWith (/) . fmap (6*) . selfDifference $ ys
--          us = fmap (2*) . zipWith (+) hs $ hs
--          vs = selfDifference bs

