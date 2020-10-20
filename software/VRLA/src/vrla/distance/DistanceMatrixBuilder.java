/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.distance;

import vrla.input.PropertiesConfiguration;

/**
 *
 * @author klin-sve
 */
public class DistanceMatrixBuilder {

    public static float[][] calculateDistanceMatrix(float rData[][]) {
        String method = PropertiesConfiguration.PROP.getProperty("method");
        if (method.equals("PEARSON")) {
            return new SquaredCorrelationDistanceMatrix().getCorrelationDistanceMatrix(rData);
        } else if (method.equals("EUCLEDIAN")) {
                return new SquaredEucledianDistanceMatrix().getEucledianDistanceMatrix(rData);
        } else if (method.equals("JACCARD")) {
                return new SquaredJaccardDistanceMatrix().getJaccardDistance(rData);
        }
        return null;
    }
}
