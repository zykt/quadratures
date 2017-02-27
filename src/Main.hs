module Main where


import Graphics.Gnuplot.Simple


data DefIntegral a = DefIntegral (a->a) a a


--instance Show DefIntegral where
--  show (DefIntegral _ start end) = "Integral from" ++ show a ++ "to" ++ show b


--alpha is power of p(x)
newton_cotes :: (Fractional a, Floating a, Num a, Enum a) => DefIntegral a -> a -> a -> a
newton_cotes (DefIntegral f start end) alpha steps =
  sum . map (\(left, right) ->
               let weight0 = ((right-start)**(1-alpha) - (left-start)**(1-alpha)) / (1-alpha)
                   weight1 = ((right-start)**(2-alpha) - (left-start)**(2-alpha)) / (2-alpha) + weight0*start
                   weigth2 = ((right-start)**(3-alpha) - (left-start)**(3-alpha)) / (3-alpha) + 2*weight1*start+ weight0*start**2
                   avg = (right - left) / 2
                   coef0 = (weigth2 - weight1*(avg + right) + weight0*avg*right) / ((avg - left)*(right - left))
                   coef1 = -(weigth2 - weight1*(left + right) + weight0*left*right) / ((avg - left)*(right - avg))
                   coef2 = (weigth2 - weight1*(avg + left) + weight0*avg*left) / ((right - left)*(right - left))
               in coef0 * f left + coef1 * f avg + coef2 * f right
            ) $ intervals start end steps

nc :: (Fractional a, Floating a, Num a, Enum a) => DefIntegral a -> a -> a
nc (DefIntegral f start end) alpha =
  let left = start
      right = end
      weight0 = ((right-start)**(1-alpha) - (left-start)**(1-alpha)) / (1-alpha)
      weight1 = ((right-start)**(2-alpha) - (left-start)**(2-alpha)) / (2-alpha) + weight0*start
      weigth2 = ((right-start)**(3-alpha) - (left-start)**(3-alpha)) / (3-alpha) + 2*weight1*start - weight0*start**2
      avg = (right - left) / 2
      coef0 = (weigth2 - weight1*(avg + right) + weight0*avg*right) / ((avg - left)*(right - left))
      coef1 = -(weigth2 - weight1*(left + right) + weight0*left*right) / ((avg - left)*(right - avg))
      coef2 = (weigth2 - weight1*(avg + left) + weight0*avg*left) / ((right - avg)*(right - left))
  in coef0 * f left + coef1 * f avg + coef2 * f right


left_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
left_sum (DefIntegral f start end) steps =
  sum . map (\(left, right) -> (right - left) * (f left)) $ intervals start end steps


average_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
average_sum (DefIntegral f start end) steps =
  sum . map (\(left, right) -> (right - left) * (f $ (left + right) / 2)) $ intervals start end steps


trapezoid_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
trapezoid_sum (DefIntegral f start end) steps =
  sum . map (\(left, right) -> (right - left) * (f right + f left) / 2) $ intervals start end steps

simpsons_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
simpsons_sum (DefIntegral f start end) steps =
  sum . map (\(left, right) -> (end - start)/steps )$ intervals start end steps

simpsons :: (Fractional a, Num a) => DefIntegral a -> a
simpsons (DefIntegral f start end) =
  (end - start) / 6 * (f start + 4 * f ((start + end) / 2) + f end)


intervals :: (Fractional a, Num a, Enum a) => a -> a -> a -> [(a, a)]
intervals start end n_intervals =
  zip intervals' $ tail intervals'
  where intervals' = [start, start + step .. end]
        step  = (end - start) / n_intervals


tests :: IO ()
tests = do
  putStrLn "hello world"
  let a = 0 :: Double
  let b = 10
  let func = \x -> x
  let integral = DefIntegral func a b
  let test_input = [2, 3, 5, 10, 100, 200]

  let test_left_sum = left_sum integral
  let test_left = map test_left_sum test_input

  let test_avg_sum = average_sum integral
  let test_right = map test_avg_sum test_input

  let test_trap_sum = trapezoid_sum integral
  let test_trap = map test_trap_sum test_input

  let test_simpson = simpsons integral

  let test_newtone_cotes_sum = newton_cotes integral (3/5)
  let test_newtone_cotes = map test_newtone_cotes_sum test_input

  putStrLn $ "steps: " ++ show test_input
  putStrLn $ "newton-cotes: " ++ show test_newtone_cotes
  putStrLn $ "left_sum: " ++ show test_left
  putStrLn $ "avg_sum: " ++ show test_right
  putStrLn $ "trap_sum: " ++ show test_trap
  putStrLn $ "simpsons: " ++ show test_simpson
  putStrLn $ "test_newton: " ++ show (nc integral (0.6))


main :: IO ()
main = do
  tests
  let func = \x -> 3.7 * cos(1.5 * x) * exp(-4*x / 3) + 2.4 * sin(4.5 * x) * exp(2*x / 3) + 4
  let a = 1.8 :: Double
  let b = 2.3
  let integral = DefIntegral func a b
  let steps = [3.0,4.0..500]

  let integral_value = 2.37880192961486

  let map_left = map (left_sum integral) steps
  let map_avg = map (average_sum integral) steps

  plotList [] $ zip steps $ map (-integral_value+) map_left
  plotList [] $ zip steps $ map (-integral_value+) map_avg
