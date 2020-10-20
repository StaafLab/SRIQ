/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.processor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.FutureTask;
import vrla.clustering.ClusterBuilder;
import vrla.input.PropertiesConfiguration;
import vrla.util.ResultTerm;

/**
 *
 * @author klin-sve
 */
public class ParallelClusteringProcess extends ProcessCreation {

    //required = distCutOffs.length
    private List<FutureTask<HashMap<String, HashMap<String, List<List<Integer>>>>>> taskList;
    private float[][] aData;

    public ParallelClusteringProcess(float[][] aData) {
        super(PropertiesConfiguration.distCutOffs.length);
        this.aData = aData;
        taskList = new ArrayList<>();
    }

    private void process(double cutOffs[]) {
        FutureTask<HashMap<String, HashMap<String, List<List<Integer>>>>> futureTask = new FutureTask<>(() -> {
            return ClusterBuilder.getQTClusters(aData, cutOffs);
        });
        taskList.add(futureTask);
        executor.execute(futureTask);
        futureTask = null;
    }

    public ResultTerm paralleProcessing() {
        for (int i = 0; i < processCount; i++) {
            double cutOffs[] = Arrays.copyOfRange(PropertiesConfiguration.distCutOffs, i == 0 ? 0 : (i * split), (i + 1) * split);
            process(cutOffs);
        }
        if (rem > 0) {
            double cutOffs[] = Arrays.copyOfRange(PropertiesConfiguration.distCutOffs, split * processCount, PropertiesConfiguration.distCutOffs.length);
            processCount += 1;
            process(cutOffs);
        }

        HashMap<String, List> ch = new HashMap<>();
        HashMap<String, List> sh = new HashMap<>();
        for (int j = 0; j < processCount; j++) {
            FutureTask<HashMap<String, HashMap<String, List<List<Integer>>>>> tasks = taskList.get(j);
            try {
                ch.putAll(tasks.get().get("ch"));
                sh.putAll(tasks.get().get("sh"));
                System.out.println("Summation Process " + (j + 1) + " results");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
        ResultTerm rt = new ResultTerm();
        rt.setCoreClusters(ch);
        rt.setSpiralClusters(sh);
        return rt;
    }

}
