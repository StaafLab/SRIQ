/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.processor;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.FutureTask;
import vrla.util.ResultTerm;

/**
 *
 * @author Sunnyveerla
 */
public class ParallelProcessing extends ProcessCreation {

    public int iterations;
    private List<FutureTask<List<ResultTerm>>> taskList;

    public ParallelProcessing(int required) {
        super(required);
        this.iterations = required;
        taskList = new ArrayList<>();
    }

    public List<ResultTerm> parallelProcess() {
        for (int i = 0; i < processCount; i++) {
            process(split);
        }
        if (rem > 0) {
            processCount += 1;
            process(rem);
        }
        List<ResultTerm> vr = new ArrayList(iterations);
        for (int j = 0; j < processCount; j++) {
            FutureTask<List<ResultTerm>> tasks = taskList.get(j);
            try {
                vr.addAll(tasks.get());
                System.out.println("Process " + (j + 1) + " results");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
        return vr;

    }

    private void process(int split) {
        FutureTask<List<ResultTerm>> futureTask = new FutureTask<>(() -> {
            return ParallelProcessing.startProcessing(split);
        });
        taskList.add(futureTask);
        executor.execute(futureTask);
        futureTask = null;
    }

    private static List<ResultTerm> startProcessing(int split) {
        List<ResultTerm> vr = new ArrayList(split);
        for (int i = 0; i < split; i++) {
            vr.add(new ExecuteClustering().getSTClusters());
            System.out.println("Iteration:\t"+i);
        }
        return vr;
    }

    /*public static String process(Parameters param, String className, String methodName) {
     try {
     ClassLoader cl = Thread.currentThread().getContextClassLoader();
     Class c = Class.forName(className, true, cl);
     Object t = c.newInstance();
     //System.out.println("In the class");
     Method[] allMethods = c.getDeclaredMethods();
     //System.out.println(allMethods.length);
     for (Method m : allMethods) {
     String mname = m.getName();
     //System.out.println(mname);
     if (mname.equals(methodName)) {
     m.setAccessible(true);
     m.invoke(t, param);
     //Object o = m.invoke(t,hParam);
     // System.out.println(mname+"\t"+m.invoke(t,hParam));
     }
     }
     } catch (Exception e) {
     e.printStackTrace();
     }
     return "Done";
     }*/
}
