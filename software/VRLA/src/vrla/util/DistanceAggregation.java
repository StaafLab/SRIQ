/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.util;

import vrla.distance.DistanceMatrixBuilder;
import vrla.input.PropertiesConfiguration;

/**
 *
 * @author klin-sve
 */
public class DistanceAggregation {

    public static int sumSquareCorrelationMatrixCount;

    public static synchronized float[][] sumSquareCorrelationMatrix(float[][] array1, float[][] array2) {
        int rlen = array1.length;
        for (int i = 0; i < rlen; i++) {
            for (int j = i + 1; j < rlen; j++) {
                array1[i][j] += array2[i][j];
                array1[j][i] = array1[i][j];
            }
        }
        return array1;
    }

    public static synchronized float[][] divideSquareCorrelationMatrixByNumber(float[][] array, int number) {
        int rlen = array.length;
        for (int i = 0; i < rlen; i++) {
            for (int j = i + 1; j < rlen; j++) {
                array[i][j] /= number;
                array[j][i] = array[i][j];
            }
        }
        return array;
    }

    public float[][] sumPermutationalResults(int permutations) {
        float[][] permutationalResults = new float[PropertiesConfiguration.TOTAL_SAMPLES][PropertiesConfiguration.TOTAL_SAMPLES];
//        System.out.println("in sum permuatations"+permutations);
        for (int i = 0; i < permutations; i++) {
            permutationalResults = sumSquareCorrelationMatrix(permutationalResults, DistanceMatrixBuilder.calculateDistanceMatrix(BaggingProcess.getDataBag()));
//            System.out.println("sumSquareCorrelationMatrixCount......................." + (++sumSquareCorrelationMatrixCount) + "\r");
        }
        return permutationalResults;
    }
}
