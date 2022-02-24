module Main where

import Interpolation.DividedDifferences (coefficients)
import Control.Monad (forM_)

xs = [-3,-2,-1, 0, 1, 2, 3, 4]
ys = [ 2, 6, 9, 2, 9, 6, 8, 7]

main :: IO ()
main = do
    putStrLn $ "Diferencias Divididas"
    putStrLn $ "valores de del eje x: " ++ show xs
    putStrLn $ "valores de del eje y: " ++ show ys
    putStrLn $ "Coeficientes obtenidos:"
    forM_ (zip [0..] $ coefficients ys xs) (\(ix,c) ->
        putStrLn $ "    c" ++ show ix ++ ": " ++ show c
        )

