module Main where

import Interpolation.CubicSpline (makeSpline,getSegments,showSegment)
import Control.Monad (forM_)

ts , ys :: Num a => [a]
ys = [ 2, 6, 9, 2, 9, 6, 8, 7]
ts = [-3,-2,-1, 0, 1, 2, 3, 4]

main :: IO ()
main = do
    putStrLn $ "Diferencias Divididas"
    putStrLn $ "valores de del eje x: " ++ show ts
    putStrLn $ "valores de del eje y: " ++ show ys
    putStrLn $ "Coeficientes obtenidos:"
    forM_ (zip [0..] . getSegments . makeSpline ys $ ts) (\(ix,seg) ->
        putStrLn $ "    S" ++ show ix ++ ": " ++ (showSegment seg "")
        )

