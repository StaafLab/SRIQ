/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.processor;

import vrla.util.ResultTerm;

/**
 *
 * @author klin-sve
 */
public class ExecuteClustering {
    public ResultTerm getSTClusters() {
        System.out.println("Getting Average Distance Matrix......");
        float cData[][] = new ParallelDistanceAggregationProcess().parallelProcess();
        System.out.println("Getting core/spiral clusters......");
        return new ParallelClusteringProcess(cData).paralleProcessing();
    }
}
