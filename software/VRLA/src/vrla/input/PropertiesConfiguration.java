/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.input;

import java.util.Arrays;
import vrla.data.FileParser;

/**
 *
 * @author klin-sve
 */
public class PropertiesConfiguration extends LoadProperties {
    public String propFilePath;
    public static String[] HEADER;
    public static int TOTAL_SAMPLES;
    public static int TOTAL_GENES;
    public static int SAMPLE_BAG_SIZE;
    public static int GENE_BAG_SIZE;
    public static int TOTAL_PERMUTATIONS;
    public static int TOTAL_ITERATIONS;
    public static int TOTAL_PROCESSORS =  Runtime.getRuntime().availableProcessors();;
    public static double[] distCutOffs;
    public static boolean SPIRAL;
    

    public PropertiesConfiguration(String propFilePath) {
        super(propFilePath);
    }

    public static void setPropertiesConfiguration() {
        FileParser fp = new FileParser(LoadProperties.PROP.getProperty("studyPath") + LoadProperties.PROP.getProperty("inFileName") + ".txt");
        HEADER = FileParser.header;
        TOTAL_SAMPLES = HEADER.length - 1;
        TOTAL_GENES = FileParser.rownames.length - 1;
        SAMPLE_BAG_SIZE = Integer.parseInt(PROP.getProperty("minClusterSize"));
        GENE_BAG_SIZE = Integer.parseInt(PROP.getProperty("minBagSize"));
        TOTAL_PERMUTATIONS = Integer.parseInt(PROP.getProperty("permutations"));
        TOTAL_ITERATIONS = Integer.parseInt(PROP.getProperty("iterations"));
        PROP.setProperty("totalSamples", ""+TOTAL_SAMPLES);
        PROP.setProperty("totalGenes", ""+ TOTAL_GENES);

        if (SAMPLE_BAG_SIZE == 0) {
            SAMPLE_BAG_SIZE = (int) Math.sqrt(TOTAL_SAMPLES);
        }

        if (GENE_BAG_SIZE == 0) {
            GENE_BAG_SIZE = (int) Math.sqrt(TOTAL_GENES);
        }

        PROP.setProperty("minClusterSize", ""+SAMPLE_BAG_SIZE);
        PROP.setProperty("minBagSize", ""+GENE_BAG_SIZE);
        
        distCutOffs = Arrays.stream(PROP.getProperty("distCutOff").trim().split(",")).mapToDouble(Double::parseDouble).toArray();
        SPIRAL = Boolean.valueOf(PROP.getProperty("spiral").trim());
    }
}
