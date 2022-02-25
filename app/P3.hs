module Main where

import Interpolation.CubicSpline (makeSpline,getSegments,showSegment)
import Data.Vector (indexed)
import Control.Monad (forM_)

ts , ys :: Num a => [a]
ts = [-3,-2,-1, 0, 1, 2, 3, 4]
ys = [ 2, 6, 9, 2, 9, 6, 8, 7]

main :: IO ()
main = do
    putStrLn $ "Spline Cúbico"
    putStrLn $ "valores de del eje x: " ++ show ts
    putStrLn $ "valores de del eje y: " ++ show ys
    putStrLn $ "Coeficientes obtenidos:"
    forM_ (indexed . getSegments $ makeSpline ts ys) (\(ix,seg) ->
        putStrLn $ "    S" ++ show ix ++ ": " ++ (showSegment seg "")
        )

