/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import java.util.ArrayList;
import vrla.input.PropertiesConfiguration;
import vrla.output.WriteResults;
import vrla.processor.ParallelProcessing;
import vrla.util.ResultTerm;

/**
 *
 * @author klin-sve
 */
public class test {

    public static void main(String[] args) {
        String propFilePath = "resources/test.properties";
        PropertiesConfiguration prop = new PropertiesConfiguration(propFilePath);
        PropertiesConfiguration.setPropertiesConfiguration();
        System.out.println(PropertiesConfiguration.TOTAL_SAMPLES);
       
//        ArrayList<ResultTerm> sr =  new ParallelProcessing(PropertiesConfiguration.TOTAL_ITERATIONS).parallelProcess();
//        WriteResults writeResults = new WriteResults(sr);
    }
}
