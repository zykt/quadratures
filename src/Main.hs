module Main where


import System.IO
import Graphics.Gnuplot.Simple
import Numeric.LinearAlgebra



data DefIntegral a = DefIntegral (a->a) a a


--instance Show DefIntegral where
--  show (DefIntegral _ start end) = "Integral from" ++ show a ++ "to" ++ show b


-- returns tuple of residue and required step
residue :: (Enum a, Field a) =>
  DefIntegral a
  -> a
  -> (DefIntegral a -> a -> a -> a)
  -> (a, a)
residue (DefIntegral f start end) alpha method =
  let quadrature step = method (DefIntegral f start end) alpha (intervals start end (floor (start-end)/step + 1))
      iter steps =
        -- find vector (J, Cm, Cm+1, ...) and take J from it
        let steps' = (3><1) $ take 3 steps
            rhs = cmap quadrature steps'
            m = find_m . toList . flatten $ rhs
            lhs = (3><1) [1, 1, 1] ||| cmap (**m) steps' ||| cmap (**(m+1)) steps'
        in (lhs <\> rhs) `atIndex` (0, 0)
      find_m [s1, s2, s3] = - log ((s3 - s2)/(s2 - s1)) / log l
      l = 2
      error = 10**(-6)
      find_residue steps =
        let res = iter steps
        in if res < error
           then (res, steps !! 2)
           else find_residue (tail steps)
  in find_residue []


cardano :: (Field a) => a -> a -> a -> a -> [a]
cardano a b c d =
  let p = c/a - b**2/(3*a**2)
      q = 2*b**3 / (27*a**3) - b*c/(3*a**2) + d/a
      phi = acos(3*q/(2*p) * sqrt(-3/p))
      y1 = 2*sqrt(-p/3) * cos(phi/3)
      y2 = 2*sqrt(-p/3) * cos(phi/3 + 2*pi/3)
      y3 = 2*sqrt(-p/3) * cos(phi/3 - 2*pi/3)
  in [y1 - b/(3*a), y2 - b/(3*a), y3 - b/(3*a)]


gauss :: (Enum a, Field a) => DefIntegral a -> a -> a -> a
gauss (DefIntegral f start end) alpha steps =
  sum . map (\(left, right) ->
               -- things here
               -- 0) substition
               let left' = end - left
                   right' = end - right
                   avg' = (right' + left') / 2
                   f' x = f $ end - x

                   -- 1) find weights
                   weight n = - (right'**(n-alpha+1) - left'**(n-alpha+1)) / (n-alpha+1)

                   --2) solve linear system for a's
                   ws_lhs = (3><3)
                     [ weight 0, weight 1, weight 2
                     , weight 1, weight 2, weight 3
                     , weight 2, weight 3, weight 4]

                   ws_rhs = (3><1)
                     [ -weight 3
                     , -weight 4
                     , -weight 5]

                   coefs_a = ws_lhs <\> ws_rhs

                   -- 3) solve 3-order equation for xs
                   [d, c, b] = toList . flatten $ coefs_a
                   xs = (1><3) $ cardano 1 b c d

                   -- 4) solve linear system for A's
                   helper_lhs = (1><3) [1, 1, 1] === xs === cmap (**2) xs
                   helper_rhs = (3><1)
                     [ weight 0
                     , weight 1
                     , weight 2]
                   coefs_A = helper_lhs <\> helper_rhs

                   -- then sum Ai*fi
                   fs = (1><3) [f' left', f' avg', f' right']
                   result_vector = fs <> coefs_A
               in result_vector `atIndex` (0, 0)
            ) $ intervals start end steps


--alpha is power of p(x)
newton_cotes :: (Enum a, Field a) => DefIntegral a -> a -> a -> a
newton_cotes (DefIntegral f start end) alpha steps =
  sum . map (\(left, right) ->
               let -- substitution
                   left' = end - left
                   right' = end - right
                   avg' = (right' + left') / 2
                   f' x = f $ end - x

                   weight n = - (right'**(n-alpha+1) - left'**(n-alpha+1)) / (n-alpha+1) 

                   weights = (3><1) [weight 0, weight 1, weight 2]

                   xs = (3><3)
                     [ 1, 1, 1
                     , left', avg', right'
                     , left'**2, avg'**2, right'**2]

                   fs = (3><1) [f' left', f' avg', f' right']

                   a_coefs = tr $ xs <\> weights

                   result_vector = a_coefs <> fs
               in result_vector `atIndex` (0, 0)
            ) $ intervals start end steps


left_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
left_sum (DefIntegral f start end) steps =
  sum . map (\(left, right) -> (right - left) * (f left)) $ intervals start end steps


average_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
average_sum (DefIntegral f start end) steps =
  sum . map (\(left, right) -> (right - left) * (f $ (left + right) / 2)) $ intervals start end steps


trapezoid_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
trapezoid_sum (DefIntegral f start end) steps =
  sum . map (\(left, right) -> (right - left) * (f right + f left) / 2) $ intervals start end steps


simpson_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
simpson_sum (DefIntegral f start end) steps =
  halfstep / 3 * (sum . map (\(left, middle, right) -> f left + 4 * f middle +  f right) $ interv)
  where step = (end - start) / steps
        halfstep = step / 2
        segment = [start + halfstep, start + 3*halfstep .. end - halfstep]
        interv = map (\x -> (x - halfstep, x, x + halfstep)) segment


simpsons :: (Fractional a, Num a) => DefIntegral a -> a
simpsons (DefIntegral f start end) =
  (end - start) / 6 * (f start + 4 * f ((start + end) / 2) + f end)


intervals :: (Fractional a, Num a, Enum a) => a -> a -> a -> [(a, a)]
intervals start end n_intervals =
  zip intervals' (tail intervals')
  where intervals' = [start, start + step .. end - step] ++ [end]
        -- add end to the end to escape problems with floating point arithmetics
        step  = (end - start) / n_intervals


tests :: IO ()
tests = do
  putStrLn "hello world"
  -- let a = 0 :: Double
  -- let b = 10
  -- let func = \x -> x
  let func = \x -> 3.7 * cos(1.5 * x) * exp(-4*x / 3) + 2.4 * sin(4.5 * x) * exp(2*x / 3) + 4
  let a = 1.8 :: Double
  let b = 2.3

  let integral = DefIntegral func a b
  let test_input = [2, 3, 5, 10, 50, 100, 200]

  let test_left_sum = left_sum integral
  let test_left = map test_left_sum test_input

  let test_avg_sum = average_sum integral
  let test_right = map test_avg_sum test_input

  let test_trap_sum = trapezoid_sum integral
  let test_trap = map test_trap_sum test_input

  let test_simpson = simpsons integral

  let test_newton_cotes_sum = newton_cotes integral (3/5)
  let test_newton_cotes = map test_newton_cotes_sum test_input

  let test_gauss_sum = gauss integral (3/5)
  let test_gauss = map test_gauss_sum test_input

  putStrLn $ "steps: " ++ show test_input
  putStrLn $ "newton-cotes: " ++ show test_newton_cotes
  -- should be 1.18515
  putStrLn $ "gauss: " ++ show test_gauss
  putStrLn $ "left_sum: " ++ show test_left
  putStrLn $ "avg_sum: " ++ show test_right
  putStrLn $ "trap_sum: " ++ show test_trap
  putStrLn $ "simpsons: " ++ show test_simpson


outputToFile :: IO()
outputToFile = do
  let func = \x -> 3.7 * cos(1.5 * x) * exp(-4*x / 3) + 2.4 * sin(4.5 * x) * exp(2*x / 3) + 4
  let a = 1.8 :: Double
  let b = 2.3
  let integral = DefIntegral func a b
  let steps = [3.0,4.0..500]

  let integral_value = 2.37880192961486

  let map_left = map (left_sum integral) steps
  let map_avg = map (average_sum integral) steps
  tests

  file <- openFile "outputnc.txt" WriteMode
  --hPutStrLn file ("integral_result = " ++ show integral_value ++ "\n")
  --hPutStrLn file ("steps = " ++ (show steps) ++ "\n")
  --hPutStrLn file ("left_sum = " ++ (show $ map (left_sum integral) steps) ++ "\n")
  --hPutStrLn file ("avg_sum = " ++ (show $ map (average_sum integral) steps) ++ "\n")
  --hPutStrLn file ("trap_sum = " ++ (show $ map (trapezoid_sum integral) steps) ++ "\n")
  hPutStrLn file ("simpson_sum = " ++ (show $ map (simpson_sum integral) steps) ++ "\n")
  --hPutStrLn file ("newton_cotes_sum = " ++ (show $ map (newton_cotes integral 0.6) steps) ++ "\n")
  hClose file


main :: IO ()
main = do
  putStrLn "This is Main.hs"
  let func = \x -> 3.7 * cos(1.5 * x) * exp(-4*x / 3) + 2.4 * sin(4.5 * x) * exp(2*x / 3) + 4
  let a = 1.8 :: Double
  let b = 2.3
  let integral = DefIntegral func a b
  let steps = [3.0,4.0..500]

  --let integral_value = 2.37880192961486
  let integral_value = 1.18515

  let map_left = map (left_sum integral) steps
  let map_avg = map (average_sum integral) steps
  let map_nc = map (newton_cotes integral 0.6) steps
  let map_gauss = map (gauss integral 0.6) steps

  --plotList [] $ zip steps $ map (-integral_value+) map_left
  --plotList [] $ zip steps $ map (-integral_value+) map_nc
  plotList [] $ zip steps $ map (-integral_value+) map_gauss
  --outputToFile
  tests
