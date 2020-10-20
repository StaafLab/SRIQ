package vrla.clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author SRIN1
 */
//Modified program and works...tested
public class SQTC {

    /**
     * @param args the command line arguments
     */
    private float fdata[][];
    private double adjustedDiameter;
    private int minimumClusterSize;
    private int number_of_items;

    public SQTC(float fdata[][], double adjustedDiameter, int minimumClusterSize) {
        this.adjustedDiameter = adjustedDiameter;
        this.minimumClusterSize = minimumClusterSize;
        this.fdata = fdata;
        this.number_of_items = fdata.length;

    }

    public synchronized List<List<Integer>> getAllClusters() {

        System.out.println("line count::::::" + number_of_items+"\t Diameter:::::"+adjustedDiameter);
        List<Integer> unassignedUniqueIDIndices = new ArrayList<>();
        for (int i = 0; i < number_of_items; i++) {
            unassignedUniqueIDIndices.add(i);
        }
//        unassignedUniqueIDIndices = Collections.synchronizedList(unassignedUniqueIDIndices);
        int fc = 0;
        int rg[] = new int[number_of_items];
        int remaining = 0;
        HashMap<String, List> hCluster = new HashMap<>();
        List<List<Integer>> allClusters = new ArrayList<>();        
        while (true) {
            List<Integer> vLargestCluster = new ArrayList<>();
            
            int largestClusterSize = 0;
            for (int q = 0; q < number_of_items; q++) {
                if (rg[q] == 1) {
                    continue;
                }
                List<Integer> mvt = new ArrayList<>();
                mvt.add(q);
                int neighbourIndex = q;
                float gd[] = new float[number_of_items];
                for (int g = 0; g < gd.length; g++) {
                    gd[g] = Float.NEGATIVE_INFINITY;
                }
                int tmp[] = new int[number_of_items];
                while (true) {
                    int indexPair = -1;
                    float bestValue = Float.POSITIVE_INFINITY;
                    //SEARCH:
                    for (int v = 0; v < number_of_items;) {
                        if (v == q || rg[v] == 1 || tmp[v] == 1) {
                            //System.out.println("out****"+v);
                            v++;
                            continue;
                        }

                        gd[v] = (float) Math.round((Math.max(fdata[v][neighbourIndex], gd[v])) * 100) / 100;

                        if (gd[v] > this.adjustedDiameter) {
                            tmp[v] = 1;
                            continue;
                            //v = 0;
                            //continue SEARCH;
                        }
                        if (gd[v] < bestValue) {
                            bestValue = gd[v];
                            indexPair = v;
                            //System.out.println(gd[v]+"\t"+d);
                            //System.out.println(indexPair);
                        }
                        v++;
                    }
                    if (indexPair == -1) {
                        break;
                    }
                    mvt.add(indexPair);
                    neighbourIndex = indexPair;
                    tmp[indexPair] = 1;
                }

                int clusterSize = mvt.size();
                if (largestClusterSize < clusterSize) {
                    largestClusterSize = clusterSize;
                    vLargestCluster.clear();
                    vLargestCluster.addAll(mvt);
                }
                //System.out.println(largestClusterSize);
            }
//            Object sp[] = vLargestCluster.toArray();
            int lcSize = vLargestCluster.size();

            if (lcSize < minimumClusterSize) {
                if (unassignedUniqueIDIndices.isEmpty()) {
                    allClusters.add(new ArrayList<Integer>());
                } else {
                    allClusters.add(unassignedUniqueIDIndices);
                }
                break;
            } else {
                List<Integer> removeItems = new ArrayList<>();
                for (int o = 0; o < lcSize; o++) {
                    int cItem = vLargestCluster.get(o);
                    rg[cItem] = 1;
//                    System.out.println("Its here");
//                    unassignedUniqueIDIndices.remove(cItem);
                    removeItems.add(cItem);
                }
                unassignedUniqueIDIndices.removeAll(removeItems);
//                 System.out.println("Its here" + unassignedUniqueIDIndices.size());
                remaining += lcSize;
                System.out.println("remaining:" + (number_of_items - remaining));
                ++fc;
                allClusters.add(vLargestCluster);
                //writeQTCluster(util, gfile, row, qtf + (fc) + "_QTCluster_" + fileName, geneListFolder + (fc) + "_QTCluster_" + fileName);
            }
            
//            System.out.println("Its here" + unassignedUniqueIDIndices.size());
        }
        return allClusters;
    }
}
