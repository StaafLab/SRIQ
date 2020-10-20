/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.processor;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.FutureTask;
import vrla.input.PropertiesConfiguration;
import vrla.util.DistanceAggregation;

/**
 *
 * @author klin-sve
 */
public class ParallelDistanceAggregationProcess extends ProcessCreation {

    private List<FutureTask<float[][]>> taskList;

    public ParallelDistanceAggregationProcess() {
        super(PropertiesConfiguration.TOTAL_PERMUTATIONS);
        taskList = new ArrayList<>();
    }

    public float[][] parallelProcess() {
        for (int i = 0; i < processCount; i++) {
            process(split);
        }
        if (rem > 0) {
            processCount += 1;
            process(rem);
        }

        float[][] permutationalResults = new float[PropertiesConfiguration.TOTAL_SAMPLES][PropertiesConfiguration.TOTAL_SAMPLES];
        for (int j = 0; j < processCount; j++) {
            FutureTask<float[][]> tasks = taskList.get(j);
            try {
                permutationalResults = DistanceAggregation.sumSquareCorrelationMatrix(permutationalResults, tasks.get());
                System.out.println("Summation Process " + (j + 1) + " results");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
       return DistanceAggregation.divideSquareCorrelationMatrixByNumber(permutationalResults, PropertiesConfiguration.TOTAL_PERMUTATIONS);
    }

    private void process(int split) {
        FutureTask<float[][]> futureTask = new FutureTask<>(() -> {
            return new DistanceAggregation().sumPermutationalResults(split);
        });
        taskList.add(futureTask);
        executor.execute(futureTask);
        futureTask = null;
    }
}
