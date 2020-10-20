/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.clustering;

import java.util.ArrayList;

/**
 *
 * @author Sunnyveerla
 */
public class QTC {

    private float fdata[][];
    private float adjustedDiameter;
    private int minimumClusterSize;
    private int number_of_items;

    public QTC(float fdata[][], float adjustedDiameter, int minimumClusterSize) {
        this.adjustedDiameter = adjustedDiameter;
        this.minimumClusterSize = minimumClusterSize;
        this.fdata = fdata;
        this.number_of_items = fdata.length;

    }

    public ArrayList<ArrayList<Integer>> getAllClusters(ArrayList<Integer> unassignedUniqueIDIndices) {
        ArrayList<ArrayList<Integer>> allClusters = new ArrayList<>();
        while (true) {

            // main work segment
            ArrayList<Integer> currentLargestCluster = getLargestCluster(unassignedUniqueIDIndices);

            allClusters.add(currentLargestCluster);

            int clusterSize = currentLargestCluster.size();

            System.out.println("# of assigned genes: " + (number_of_items - unassignedUniqueIDIndices.size()));

            System.out.println("# of genes not yet assigned: " + unassignedUniqueIDIndices.size());

            System.out.println("# of clusters formed: " + allClusters.size());

            System.out.println("size of last cluster formed: " + currentLargestCluster.size());

            unassignedUniqueIDIndices.removeAll(currentLargestCluster);

            // terminating conditions
            if (clusterSize < minimumClusterSize) { // throw the leftovers into the last cluster.

                currentLargestCluster.addAll(unassignedUniqueIDIndices);

                return allClusters;

            } else if (unassignedUniqueIDIndices.size() == 0) { // insert an empty cluster to indicate it came out evenly.

                allClusters.add(new ArrayList());

                return allClusters;

            }

        } // end while true

    } // end getting clusters

    private ArrayList<Integer> getLargestCluster(ArrayList<Integer> unassignedUIDIndices) {

        ArrayList<Integer> currentCluster;

        ArrayList<Integer> largestClusterTies = new ArrayList<>();

        int largestClusterSize = 0;

        for (int i = 0; i < unassignedUIDIndices.size(); i++) {

            ArrayList<Integer> tempUnassigned = (ArrayList<Integer>) unassignedUIDIndices.clone(); // need to clone because getClusterForAGene mangles the unassigned indices array.

            Integer seedCandidate = tempUnassigned.remove(i);

            currentCluster = getClusterForAGene(seedCandidate, tempUnassigned);

           /* if (currentCluster.size() == largestClusterSize) { // if a tie

                largestClusterTies.add(currentCluster);

            } else*/
            if (currentCluster.size() > largestClusterSize) {

                largestClusterTies.clear();

                largestClusterTies.addAll(currentCluster);

                largestClusterSize = currentCluster.size();

            } // if we have a new record

        } // for each possible seed gene

        //int randCluster = (int) (Math.random() * largestClusterTies.size());

        return largestClusterTies;

    }

    /**
     *
     * This returns the cluster that a gene would produce. This has specific
     * permission to
     *
     * mangle the vector of indices passed to it.
     *
     *
     *
     * @param candidateIndices does NOT contain an entry for seedIndex.
     *
     */
    private ArrayList<Integer> getClusterForAGene(Integer seedIndex, ArrayList<Integer> candidateIndices) {

        ArrayList<Integer> cluster = new ArrayList();

        cluster.add(seedIndex);

        Integer MostRecentAdditionI = seedIndex;

        int mostRecentAdditioni = seedIndex.intValue();

        // the potential diameter for the cluster if any one gene is added. The worst of the distances.
        float[] geneDiameterSoFar = new float[number_of_items]; // indexed by absolute indices so it can be non-dynamic and primitive

        for (int local = 0; local < candidateIndices.size(); local++) {
            geneDiameterSoFar[(candidateIndices.get(local)).intValue()] = Float.NEGATIVE_INFINITY;
        }

        while (true) { // exit condition near bottom

            int bestLocalIndex = -1;

            float bestDistance = Float.POSITIVE_INFINITY; // best of the worst distances.

            //    initialize the array only for those genes that we're using.
            // main loop
            CANDIDATE_SEARCH:
            for (int local = 0; local < candidateIndices.size();) { // increment local only if not deleting!

                int i = (candidateIndices.get(local)).intValue(); // the absolute index at the local index

                // compare this gene to each gene already in the cluster. Keep track of worst match.
                geneDiameterSoFar[i] = Math.max(fdata[i][mostRecentAdditioni], geneDiameterSoFar[i]); // potential diameter
                
                
                if (geneDiameterSoFar[i] > adjustedDiameter) { // never check this gene again

                    candidateIndices.remove(local); // local now points to the element after the one it used to.

                    continue CANDIDATE_SEARCH;

                } // if the candidate gene is disqualified now

                if (geneDiameterSoFar[i] < bestDistance) { // if this worst match is the best one so far, make it the leader.

                    bestDistance = geneDiameterSoFar[i];

                    bestLocalIndex = local;

                }

                local++;

            } // for each candidate gene

            if (bestLocalIndex == -1) {
                break; // if no candidates can join cluster, stop adding them already!
            }

            MostRecentAdditionI = candidateIndices.remove(bestLocalIndex); // otherwise, add the best candidate.

            mostRecentAdditioni = MostRecentAdditionI.intValue();

            cluster.add(MostRecentAdditionI);

        } // while can add genes

        return cluster;

    }

}
