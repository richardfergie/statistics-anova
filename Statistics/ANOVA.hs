module Statistics.ANOVA (anova, fisherLSD) where
import Data.List
import Statistics.Distribution hiding (mean, variance)
import Statistics.Distribution.FDistribution
import Statistics.Distribution.StudentT
import Statistics.Sample
import qualified Data.Vector as V

{-# INLINE treatments #-}
treatments :: (Ord a,Eq a) => V.Vector (a, Double) -> V.Vector a
treatments input = V.fromList $! nub $! V.toList $! V.map fst input --not good

treatmentMeans :: (Ord a, Eq a) => V.Vector a -> V.Vector (a, Double) -> V.Vector (a,Double) 
treatmentMeans treatments input = V.map (\x -> x input) meanFuncs 
  where meanFuncs = V.map treatmentMean $! treatments 

treatmentMean :: (Ord a, Eq a) => a -> V.Vector (a,Double) -> (a,Double)
treatmentMean treatment input = (treatment, avg) 
  where avg = mean $! V.map snd $! filterTreatment treatment $! input
        
treatmentSize treatment input = V.length $! filterTreatment treatment $! input        
        
{-# INLINE filterTreatment #-}
filterTreatment :: (Ord a, Eq a) => a -> V.Vector (a, Double) -> V.Vector (a, Double)
filterTreatment treatment = V.filter (\(x,y) -> x == treatment)

sampleMean input = mean $! V.map snd input        

ssTreatment :: (Ord a, Eq a) => V.Vector a -> V.Vector (a, Double) -> Double
ssTreatment treatments input = V.sum $! V.zipWith (\size mean -> size*(mean-sMean)*(mean-sMean)) tSizes tMeans  
  where sMean = sampleMean input
        tMeans = V.map snd $! treatmentMeans treatments input
        tSizes = V.map (\x -> fromIntegral $! x input) $! V.map treatmentSize $! treatments --input

ssError :: (Ord a, Eq a) => V.Vector (a, Double) -> Double
ssError input = V.sum $! V.map (\x -> x input) $! V.map treatmentError $! treatments input

treatmentError treatment input = V.sum $! V.map (\x -> (x-tMean)*(x-tMean)) filteredInput
  where filteredInput = V.map snd $! filterTreatment treatment input 
        tMean = snd $! treatmentMean treatment input

fValue :: (Ord a, Eq a) => V.Vector a -> V.Vector (a, Double) -> Double
fValue treatments input = ((ssTreatment treatments input)/(ssError input)) * ((sampleSize - numberOfTreatments)/(numberOfTreatments -1))
  where numberOfTreatments = fromIntegral $! V.length $! treatments 
        sampleSize = fromIntegral $! V.length input
        
fDist treatments input = fDistribution (numberOfTreatments-1) (sampleSize-numberOfTreatments)      
  where numberOfTreatments = fromIntegral $! V.length treatments 
        sampleSize = fromIntegral $! V.length input
        
pValue input = 1 - cumulative (fDist t input) (fValue t input)
  where t = treatments input

anova pvalue input
  | pValue input < pvalue = True
  | otherwise = False                          

tTest :: (Ord a, Eq a) => Double -> a -> a -> V.Vector (a, Double) -> Bool
tTest pvalue treatment1 treatment2 inputs = abs (tTestT treatment1 treatment2 inputs) > critical
  where critical = quantile (studentT (degreesOfFreedom treatment1 treatment2 inputs)) (1 - pvalue/2)

tTestT :: (Ord a, Eq a) => a -> a -> V.Vector (a, Double) -> Double
tTestT treatment1 treatment2 inputs = (mean1-mean2)/sqrt(var1/size1+var2/size2)
  where mean1 = snd $! treatmentMean treatment1 inputs
        mean2 = snd $! treatmentMean treatment2 inputs
        size1 = fromIntegral $! V.length $! filterTreatment treatment1 inputs
        size2 = fromIntegral $! V.length $! filterTreatment treatment2 inputs
        var1 = variance $! V.map snd $! filterTreatment treatment1 inputs
        var2 = variance $! V.map snd $! filterTreatment treatment2 inputs
        
degreesOfFreedom
  :: (Ord a, Eq a) => a -> a -> V.Vector (a, Double) -> Double
degreesOfFreedom treatment1 treatment2 inputs = numerator / denominator
  where size1 = fromIntegral $! V.length $! filterTreatment treatment1 inputs
        size2 = fromIntegral $! V.length $! filterTreatment treatment2 inputs
        var1 = variance $! V.map snd $! filterTreatment treatment1 inputs
        var2 = variance $! V.map snd $! filterTreatment treatment2 inputs
        numerator = (var1/size1+var2/size2)^2
        denominator = (((var1/size1)^2)/(size1-1)) + (((var2/size2)^2)/(size2-1))

noDifferenceInMean pvalue treatment input = (:) treatment $!
                                            map fst $!
                                            filter (\(x,y) -> not y) $! 
                                            map (\x -> (x, tTest pvalue treatment x input)) othertreatments
  where othertreatments = delete treatment $! V.toList $! treatments input
        
equivalenceList pvalue input = nub $! V.toList $! V.map (\x -> sort $! noDifferenceInMean pvalue x input) $! treatments input

fisherLSD :: (Eq a, Ord a) => Double -> V.Vector (a, Double) -> [[a]]
fisherLSD pvalue input = case (anova pvalue input) of
  False -> [V.toList $! treatments input]
  True -> equivalenceList pvalue input
