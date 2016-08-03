import unittest
import numpy as np
import startingConformationsCalculator
import clustering

class TestStartingConformationsCalculator(unittest.TestCase):
    def testDivideTrajAccordingToWeights(self):
        startingConfCalculator = startingConformationsCalculator.StartingConformationsCalculator() 
        weights1 = [0.5, 0.2, 0.2, 0.1]
        trajToDistribute1 = 12
        degeneracy1 = startingConfCalculator.divideTrajAccordingToWeights(weights1, trajToDistribute1)

        golden1 = [6, 2, 3, 1]
        self.assertEqual(degeneracy1, golden1)

        falsegolden1 = [6, 2, 3, 2]
        self.assertNotEqual(degeneracy1, falsegolden1)


    def testDivideInverselyProportionalToArray(self):
        startingConfCalculator = startingConformationsCalculator.StartingConformationsCalculator() 
        weights2 = np.array([0.5, 0.2, 0.2, 0.1])
        trajToDistribute2 = 12
        degeneracy2 = startingConfCalculator.divideInverselyProportionalToArray(weights2, trajToDistribute2)

        golden2 = [1, 3, 3, 5]

        test2Passed = True
        self.assertEqual(degeneracy2, golden2)

        falsegolden2 = [1, 3, 3, 6]
        self.assertNotEqual(degeneracy2, falsegolden2)

    def testInverselyProportionalToPopulationCalculator(self):
        #test 3
        clusters = clustering.Clusters()
        sizes = [6,2,3,1]
        for size in sizes:
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            clusters.addCluster(cluster)

        inverselyProp = startingConformationsCalculator.InverselyProportionalToPopulationCalculator()
        clusteringParams = None
        trajs = 10
        degeneracy3 = inverselyProp.calculate(clusters.clusters, trajs, clusteringParams)
        golden3 = [1,2,2,5]

        self.assertEqual(degeneracy3, golden3)

        falsegolden3 = [1,2,2,6]
        self.assertNotEqual(degeneracy3, falsegolden3)


    def testEpsilonCalculator(self):
        epsilon = startingConformationsCalculator.EpsilonDegeneracyCalculator()
        params = startingConformationsCalculator.SpawningParams()
        params.epsilon = 0.5

        clusters = clustering.Clusters()
        sizes = [6,2,3,1]
        energies = [-4,-2,-2,-1]
        for size, energy in zip(sizes, energies):
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            cluster.metric = energy
            clusters.addCluster(cluster)

        trajs = 20
        degeneracy4 = epsilon.calculate(clusters.clusters, trajs, params)
        golden4 = np.array([7,4,4,5])
        np.testing.assert_array_equal(degeneracy4, golden4)

        #test 5
    def testSameWeightDegeneracyCalcuulator(self):
        sameWeightDegCalculator = startingConformationsCalculator.SameWeightDegeneracyCalculator()
        params = None

        clusters = clustering.Clusters()
        sizes = [6,2,3,1]
        for size in sizes:
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            clusters.addCluster(cluster)

        trajs = 10
        degeneracy5 = sameWeightDegCalculator.calculate(clusters.clusters, trajs, params)
        golden5 = [1, 1, 1, 1]

        self.assertEqual(degeneracy5, golden5)

        falseGolden5 = [1, 1, 1, 2]
        self.assertNotEqual(degeneracy5, falseGolden5)

def main():
    return unittest.main(exit=False)

if __name__ == '__main__':
    main()