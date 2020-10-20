/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.processor;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import vrla.input.PropertiesConfiguration;

/**
 *
 * @author klin-sve
 */
public class ProcessCreation {

    public int processCount = PropertiesConfiguration.TOTAL_PROCESSORS;
    public ExecutorService executor;
    public final int split;
    public final int rem;

    public ProcessCreation(int required) {
        if (required <= processCount) {
            processCount = required;
        }
        split = required / processCount;
        rem = required % processCount;

        executor = Executors.newFixedThreadPool(processCount);
    }

}
