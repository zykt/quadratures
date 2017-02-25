module Main where


data DefIntegral a = DefIntegral (a->a) a a

--instance Show DefIntegral where
--  show (DefIntegral _ start end) = "Integral from" ++ show a ++ "to" ++ show b

general_sum :: (Fractional a, Num a, Enum a) => ((a, a) -> a) -> DefIntegral a -> a ->  a
general_sum (DefIntegral f start end) steps sum_func =
  sum . map sum_func $ intervals start end steps

left_sum' :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
left_sum' = general_sum (\(left, right) -> (right - left) * (f left))

left_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
left_sum (DefIntegral f start end) steps =
  integral_sum
  where integral_sum = sum . map (\(left, right) -> (right - left) * (f left)) $ intervals start end steps

average_sum :: (Fractional a, Num a, Enum a) => DefIntegral a -> a -> a
average_sum (DefIntegral f start end) steps =
  integral_sum
  where integral_sum = sum . map (\(left, right) -> (right - left) * (f $ (left + right) / 2)) $ intervals start end steps


uniform_steps :: (Fractional a, Num a, Enum a) => a -> a -> a -> [a]
uniform_steps start end steps =
  [start, start + step .. end]
  where step = (end - start) / steps

intervals :: (Fractional a, Num a, Enum a) => a -> a -> a -> [(a, a)]
intervals start end n_intervals =
  zip intervals' $ tail intervals'
  where intervals' = [start, start + step .. end]
        step  = (end - start) / n_intervals

main :: IO ()
main = do
  putStrLn "hello world"
  let a = 0 :: Double
  let b = 10
  --let steps = 100
  --let integral1 = left_sum (DefIntegral (\x->x) a b) steps
  let func = \x -> x*x
  let integral = DefIntegral func a b
  let test_input = [2, 3, 5, 10, 100, 200]
  let test_left_sum = left_sum integral
  let test_left = map test_left_sum test_input
  let test_left'_sum = left_sum' integral
  let test_left' = map test_left'_sum test_input
  let test_avg_sum = average_sum integral
  let test_right = map test_avg_sum test_input
  --putStrLn $ show $ intervals a (b - (b-a)/steps) steps
  --putStrLn $ show integral1
  putStrLn $ "steps: " ++ show test_input
  putStrLn $ "left_sum: " ++ show test_left
  putStrLn $ "avg_sum: " ++ show test_right
  
