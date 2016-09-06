import os
import trajectories
import helper
import msm
import tpt
import json
import matplotlib.pyplot as plt

class PrepareMSM:
    def __init__(self,numClusters,trajectoryFolder, trajectoryBasename ):

        self.discretizedFolder = "discretized"
        self.clusterCentersFile = os.path.join(self.discretizedFolder, "clusterCenters.dat")
        self.discTraj = os.path.join(self.discretizedFolder, "%s.disctraj")
        self.clusteringFile = "clustering_object.pkl"
        if os.path.exists(self.clusteringFile):
            self.cl = helper.loadClustering(self.clusteringFile)
        else:
            self.loadTrajectories(numClusters, trajectoryFolder,
                                  trajectoryBasename)

    def loadTrajectories(self, numClusters, trajectoryFolder,
                         trajectoryBasename):
        print "Loading trajectories..."
        self.x = trajectories.loadCOMFiles(trajectoryFolder, trajectoryBasename)

        #cluster & assign
        print "Clustering data..."
        self.cl = trajectories.clusterTrajectories(self.x, numClusters)
        #write output
        print "Writing clustering data..."
        helper.makeFolder(self.discretizedFolder)
        helper.writeClusterCenters(self.cl, self.clusterCentersFile)
        print "Saving clustering object..."
        helper.saveClustering(self.cl, self.clusteringFile)

    def getClusteringObject(self):
        return self.cl

class MSM:
    def __init__(self, cl, lagtimes, numPCCA, itsOutput=None, numberOfITS=-1,
                 itsErrors=None, error_estimationCK=None):
        self.MSMFile = "MSM_object.pkl"
        if os.path.exists(self.MSMFile):
            self.MSM_object = helper.loadMSM(self.MSMFile)
        else:
            lagtime = self.calculateITS(cl,lagtimes,itsOutput, numberOfITS, itsErrors)
            self.createMSM(cl, lagtime)
            self.check_connectivity()
            self.PCCA(numPCCA)
            print "Saving MSM object..."
            helper.saveMSM(self.MSM_object)
            self.performCPTest(error_estimationCK)
    
    def getMSM_object(self):
        return self.MSM_object

    def performCPTest(self,error_estimationCK=None):
        #Chapman-Kolgomorov validation
        nsetsCK = len(self.MSM_object.metastable_sets)
        print ("Performing Chapman-Kolmogorov validation with the %d sets from the "
            "PCCA, when it's done will prompt for the validity of the model...")%nsetsCK
        membershipsCK = self.MSM_object.metastable_memberships
        CKObject = msm.ChapmanKolmogorovTest(self.MSM_object,
                                            nsetsCK,memberships=membershipsCK,
                                            error_estimation=error_estimationCK)
        msm.plotChapmanKolmogorovTest(CKObject)
        plt.show()

    def PCCA(self, numPCCA):
        #PCCA
        print "Calculating PCCA cluster with %d sets..."%numPCCA
        self.MSM_object = msm.calculatePCCA(self.MSM_object, numPCCA)

    def check_connectivity(self):
        #connectivity
        print "Checking connectivity of the MSM..."
        if msm.is_connected(self.MSM_object):
            print "The MSM estimated is fully connected"
        else:
            print "The MSM estimated is not fully connected"
            unconnected_sets = self.MSM_object.connected_sets
            print "The MSM estimated has %d connected sets with sizes:" % len(unconnected_sets)
            for index, uncon_set in enumerate(unconnected_sets):
                print "Set %d has %d elements" % (index, uncon_set.size)


    def createMSM(self, cl, lagtime):
        #estimation
        print "Estimating MSM with lagtime %d..."%lagtime
        self.MSM_object = msm.estimateMSM(cl.dtrajs, lagtime)

    def calculateITS(self, cl, lagtimes,itsOutput=None, numberOfITS=-1,
                     itsErrors=None):
        is_converged = False
        #its
        print ("Calculating implied time-scales, when it's done will prompt for "
                "confirmation on the validity of the lagtimes...")
        while not is_converged:
            its_object = msm.calculateITS(cl.dtrajs, lagtimes, itsErrors)
            plot_its = msm.plotITS(its_object, itsOutput, numberOfITS)
            plt.show()
            while True:
                convergence_answer = raw_input("Has the ITS plot converged?[y/n] ")
                convergence_answer.rstrip()
                convergence_answer = convergence_answer or "y" #Making yes the default
                #answer
                if convergence_answer.lower() == "y" or convergence_answer.lower() == "yes":
                    is_converged = True
                    lagtime_str = raw_input("Please input the lagtime to construct the MSM: ")
                    lagtime = int(lagtime_str.rstrip())
                    break
                elif convergence_answer.lower() == "n" or convergence_answer.lower() == "no":
                    break
                else:
                    print "Answer not valid. Please answer yes or no"
            if not is_converged:
                new_lagtimes = raw_input("Do you want to define new lagtimes or add to the previous?[add(a)/new(n)] ")
                new_lagtimes.rstrip()
                if new_lagtimes.lower() == "add" or new_lagtimes.lower() == "a":
                    lag_list = raw_input("Please input the lagtimes you want to add separated by a space: ")
                    lag_list.rstrip()
                    lagtimes.extend(map(int,lag_list.split(" ")))
                elif new_lagtimes.lower() == "new" or new_lagtimes.lower() == "n":
                    lag_list = raw_input("Please input the new lagtimes separated by a space: ")
                    lag_list.rstrip()
                    lagtimes = map(int,lag_list.split(" "))
                lagtimes.sort()
        return lagtime

class TPT:
    def __init__(self, MSM_object, cl, outfile_fluxTPT, state_labels):
        #Identify relevant sets for TPT
        print ("Plotting PCCA sets, identify the sets that will serve as source and sink "
            "for TPT...")
        msm.plot_PCCA_clusters(cl, MSM_object)
        plt.show()
        SetA_index = int(raw_input("Please input index of Set A(source):"))
        SetB_index = int(raw_input("Please input index of Set B(sink):"))
        SetA, SetB = tpt.selectTPTSets(MSM_object, SetA_index, SetB_index)
        print "Creating TPT object..."
        self.TPT_object = tpt.createTPT(MSM_object, SetA, SetB)
        print "Coarsing TPT for visualization..."
        self.coarseTPT_object = tpt.coarseTPT(self.TPT_object, MSM_object)
        print "Plotting TPT flux diagram..."
        flux_figure = tpt.plotTPT(self.coarseTPT_object, state_labels=state_labels,
                                outfile=outfile_fluxTPT)
        plt.show()
        print "Writing the main properties of the TPT in the tpt/ folder..."
        tpt.writeTPTOutput(self.coarseTPT_object)

    def getTPTObject(self):
        return self.TPT_object

    def getCoarseTPTObject(self):
        return self.coarseTPT_object

def readParams(control_file):
    """Reads a JSON control file for the paramaters.
    List of parameters:
        trajectoryFolder : string, folder where the trajectory data is
        trajectoryBasename : string, names of the trajectory files is
        numClusters : int, number of clusters to partition the data
        lagtimes: list of ints, lagtimes to calculate the implied time-scales
        numPCCA : int, number of sets to perform PCCA
        itsOutput : string (optional), name of the file where to store the implied time-scales plot
        numberOfITS : int (optional), number of eigenvalues to include in the ITS plot
        itsErrors : string (optional), 'bayes' to include error estimation in ITS plot
        error_estimationCK : bool (optional), wether to include error estimation in CK test
        state_labels : list of strings (optional), list of labels for the sets in the coarse TPT flux plot
        outfile_fluxTPT : string, name of the file to store the coarse TPT flux plots
        """
    with open(control_file, "r") as f:
        paramsJSON = json.load(f)
    return paramsJSON
