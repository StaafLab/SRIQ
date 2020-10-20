/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.output;

import Plots.Plots;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import vrla.clustering.ClusterBuilder;
import vrla.clustering.STClusterAnalysis;
import vrla.input.PropertiesConfiguration;
import vrla.util.ResultTerm;
import vrla.util.ToolBox;

/**
 *
 * @author Sunnyveerla
 */
public class WriteResults {

    private boolean spiral;
    private String tagName;
    private String header[];
    private double d;
    private int minSamples;
    private float avgTotalSamplesInvolved;
    private List<ResultTerm> vr;
    ToolBox util = new ToolBox();
    HashMap<String, List<List<Integer>>> hFinalClusters = new HashMap();
    private String keyValue;

    public WriteResults(List<ResultTerm> vr) {
        this.spiral = PropertiesConfiguration.SPIRAL;
        this.vr = vr;
        String outPath = PropertiesConfiguration.PROP.getProperty("outPath");
        String qcOutFileName = PropertiesConfiguration.TOTAL_PERMUTATIONS + "_QC";
        String qualityAnalysisFolder = outPath + PropertiesConfiguration.TOTAL_PERMUTATIONS + "/QC_Spiral(" + false + ")/";
        writeClusters(qualityAnalysisFolder, false);
        if (spiral) {
            qualityAnalysisFolder = outPath + PropertiesConfiguration.TOTAL_PERMUTATIONS + "/QC_Spiral(" + PropertiesConfiguration.SPIRAL + ")/";
            writeClusters(qualityAnalysisFolder, spiral);
        }
    }

    private void writeClusters(String qualityAnalysisFolder, boolean spiral) {
        HashMap<String, HashMap<String, Integer>> clustFreq = new HashMap();
        HashMap<String, HashMap<String, List<List<Integer>>>> clustInfo = new HashMap();
        int maxClusters = 0;
        for (int i = 0; i < PropertiesConfiguration.TOTAL_ITERATIONS; i++) {
            ResultTerm rt = vr.get(i);
            HashMap<String, List> hr = rt.getClusters(spiral);
            if (hr.size() == 0) {
                break;
            }
            for (double d : PropertiesConfiguration.distCutOffs) {
                //Frequency of clusters
                if (!hr.containsKey("" + d)) {
                    continue;
                }

                int num_of_clusters = (hr.get("" + d).size() - 1);

                //System.out.println(num_of_clusters);
                if (maxClusters < num_of_clusters) {
                    maxClusters = num_of_clusters;
                    //System.out.println("maxClusters**** "+(i+1));
                }
                List<List<Integer>> vctmp = hr.get("" + d);
                //remove the unassigned clusters                
                vctmp.remove(vctmp.get(num_of_clusters));

                if (!clustFreq.containsKey("" + d)) {
                    HashMap<String, Integer> tmp = new HashMap();
                    tmp.put("" + num_of_clusters, 1);
                    clustFreq.put("" + d, tmp);
                    HashMap<String, List<List<Integer>>> ctmp = new HashMap();
                    ctmp.put("" + num_of_clusters, vctmp);
                    clustInfo.put("" + d, ctmp);
                } else {
                    HashMap<String, Integer> tmp = clustFreq.get("" + d);
                    HashMap<String, List<List<Integer>>> ctmp = clustInfo.get("" + d);
                    if (!tmp.containsKey("" + num_of_clusters)) {
                        tmp.put("" + num_of_clusters, 1);
                        clustFreq.put("" + d, tmp);
                        ctmp.put("" + num_of_clusters, vctmp);
                        clustInfo.put("" + d, ctmp);
                    } else {
                        tmp.put("" + num_of_clusters, tmp.get("" + num_of_clusters) + 1);
                        clustFreq.put("" + d, tmp);
                        vctmp.addAll(ctmp.get("" + num_of_clusters));
                        ctmp.put("" + num_of_clusters, vctmp);
                        clustInfo.put("" + d, ctmp);
                    }
                }
            }
        }

//QC results
        new File(qualityAnalysisFolder).mkdirs();
        //Write Cluster frequencies results
        String[] results = new String[maxClusters + 1];
        String[] perresults = new String[maxClusters + 1];
        results[0] = "Num_Of_Clusters";
        perresults[0] = "Num_Of_Clusters";
        for (int i = 1; i <= maxClusters; i++) {
            results[i] = "" + i;
            perresults[i] = "" + i;
        }

        Object keys[] = clustFreq.keySet().toArray();
        for (Object k : keys) {
            results[0] += "\t" + k;
            perresults[0] += "\t" + k;
            HashMap<String, Integer> tmp = clustFreq.get(k);
            int cs = tmp.size();
            for (int c = 1; c <= maxClusters; c++) {
                this.keyValue = "" + k + "_" + c;
                if (tmp.containsKey("" + c)) {
                    String dPath = qualityAnalysisFolder + "dist(" + k + ")/DistanceMatrices/";
                    new File(dPath).mkdirs();
                    float[][] cData = getCoocurranceSummary(clustInfo.get(k).get("" + c));
                    //writeDistanceMatrix(cData, dPath + "C" + c + "_non_normalized_distanceMatrix.txt");                    
                    float[][] sData = util.substract(util.normalize(cData), (float) 1.0);
                    writeDistanceMatrix(sData, dPath + "C" + c + "_normalized_distanceMatrix.txt");
                    String fPath = qualityAnalysisFolder + "dist(" + k + ")/" + c + "/";
                    new File(fPath).mkdirs();
                    this.tagName = "" + c + "C";
                    List<List<Integer>> finalVC = new ArrayList();
                    prepareForWritingClusters(finalVC = getFinalClusters(sData, c, this.minSamples), sData, fPath, spiral);
                    if (!spiral) {
                        hFinalClusters.put(this.keyValue, finalVC);
                    }
                    results[c] += "\t" + tmp.get("" + c).intValue() + "_SB_" + this.d + "_PS_" + util.roundDouble(this.avgTotalSamplesInvolved / PropertiesConfiguration.TOTAL_SAMPLES, 2);
                    perresults[c] += "\t" + tmp.get("" + c).intValue() + "(" + (util.roundDouble(this.avgTotalSamplesInvolved / PropertiesConfiguration.TOTAL_SAMPLES, 2)) + ")";
                    cData = null;
                    sData = null;
                } else {
                    results[c] += "\t" + 0;
                    perresults[c] += "\t" + 0;
                }
            }
        }
        String clusterFreqFile = qualityAnalysisFolder + PropertiesConfiguration.PROP.getProperty("studyName") + "_" + PropertiesConfiguration.TOTAL_PERMUTATIONS + "_Clusters_Frequencies.txt";
        util.writeStringArrayToFile(results, clusterFreqFile);
        new Plots().makeConsistencyPlot(clusterFreqFile, PropertiesConfiguration.PROP.getProperty("studyName")+"_"+PropertiesConfiguration.TOTAL_PERMUTATIONS+"Per_"+PropertiesConfiguration.GENE_BAG_SIZE+"_genes");
    }

    private void prepareForWritingClusters(List<List<Integer>> clusterVector, float[][] fdata, String outPath, boolean spiral) {
        int row[] = new int[1];
        int num_of_clusters = clusterVector.size() - 1;//-1 because the last vector contains unassigned items        
        if (spiral) {
            if (num_of_clusters >= 2) {
                clusterVector = ClusterBuilder.getSpiralClusters(clusterVector, fdata);
            }
        }
        for (int i = 0; i < clusterVector.size(); i++) {
            int sp[] = clusterVector.get(i).stream().mapToInt(Integer::new).toArray();
            row = new int[sp.length + 1];
            row[0] = 0;

            for (int j = 1; j < row.length; j++) {
                row[j] = sp[j - 1] + 1;
            }
            if (i == num_of_clusters) {
                String fpath = new File(outPath).getParent() + "/Unassigned/";
                new File(fpath).mkdir();
                writeSTCluster(row, fpath + tagName + "_Unassigned");
            } else {
                int cn = i + 1;
                if (spiral) {
                    cn = getSpiralClusterNumbers(sp);
                    //System.out.println((i + 1) + "replaced by " + cn);
                }
                writeSTCluster(row, outPath + cn + "_" + tagName);
            }
        }
    }

    private int getSpiralClusterNumbers(int co[]) {
        int cn = 0;
        int m = 0;
        List<List<Integer>> coreCluster = hFinalClusters.get(keyValue);
        for (int i = 0; i < coreCluster.size() - 1; i++) {
            int sp[] = coreCluster.get(i).stream().mapToInt(Integer::new).toArray();
            List<Integer> asp = new ArrayList(Arrays.asList(sp));
            List<Integer> aco = new ArrayList(Arrays.asList(co));
            aco.retainAll(asp);
            int s = aco.size();
            if (s > m) {
                m = s;
                cn = i + 1;
            }
        }
        return cn;
    }

    private void writeSTCluster(int[] row, String geneListFolder) {
        String clust[] = util.getDataFromStringArray(PropertiesConfiguration.HEADER, row);
        try {            
            BufferedWriter bwg = new BufferedWriter(new FileWriter(geneListFolder + "_Samples.txt"));
            for (int i = 0; i < clust.length; i++) {
                String sp[] = clust[i].split("\t");
                bwg.write(sp[0] + "\n");
                bwg.flush();
            }
            bwg.flush();
            bwg.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void writeDistanceMatrix(float[][] results, String outFile) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
            String items[] = PropertiesConfiguration.HEADER;
            for (String item : items) {
                bw.write(item + "\t");
            }
            bw.write("\n");
            bw.flush();
            int totalItems = items.length - 1;
            System.out.println(totalItems);
            for (int i = 0; i < totalItems; i++) {
                bw.write(items[i + 1] + "\t");
                for (int j = 0; j < totalItems; j++) {
                    bw.write(results[i][j] + "\t");
                }
                bw.write("\n");
                bw.flush();
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private float[][] getCoocurranceSummary(List<List<Integer>> clusters) {
        int cSize = clusters.size();
        int ms = Integer.MAX_VALUE;
        int total = 0;
        float xData[][] = new float[PropertiesConfiguration.TOTAL_SAMPLES][PropertiesConfiguration.TOTAL_SAMPLES];
        for (int c = 0; c < cSize; c++) {
            List<Integer> vc = clusters.get(c);            
            int clen = vc.size();
            total += clen;
            if (clen < ms) {
                ms = clen;
            }
            for (int i = 0; i < clen; i++) {
                for (int j = i + 1; j < clen; j++) {
                    int x = vc.get(i);
                    int y = vc.get(j);                    
                    xData[x][y] += 2;
                    xData[y][x] = xData[x][y];
                }
            }
        }
        this.minSamples = ms;
        this.avgTotalSamplesInvolved = ((float) total / PropertiesConfiguration.TOTAL_ITERATIONS);
        System.out.println(this.minSamples);
        return xData;
    }

    private List<List<Integer>> getFinalClusters(float[][] sData, int num_of_clusters, int cSize) {
        System.out.println("Get Final Clusters....");        
        double d = 0.0;
        int percentiles = 20;
        int num_of_items = sData.length;
        List vc = new ArrayList(1);
        for (int p = 0; p < percentiles; p++) {
            List<List<Integer>> clusterVector = new ArrayList<>();
            clusterVector = new STClusterAnalysis(sData, d, cSize).getClusterAnalysisClassObject().getAllClusters();
            if (num_of_clusters == (clusterVector.size() - 1)) {
                vc = clusterVector;
                break;
            }
            d = util.roundDouble((d + 0.05), 2);
        }
        this.d = util.roundDouble((1 - d), 2);
        return vc;
    }
}