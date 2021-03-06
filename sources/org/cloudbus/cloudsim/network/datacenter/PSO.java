package org.cloudbus.cloudsim.network.datacenter;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sourceforge.jswarm_pso.SwarmGSA_GA;
import net.sourceforge.jswarm_pso.Myswarm2;
import net.sourceforge.jswarm_pso.SwarmGSA;
import net.sourceforge.jswarm_pso.Neighborhood;
import net.sourceforge.jswarm_pso.Neighborhood1D;
import net.sourceforge.jswarm_pso.Swarm;
import net.sourceforge.jswarm_pso.example_2.SwarmShow2D;
import org.cloudbus.cloudsim.Datacenter;
import org.cloudbus.cloudsim.network.datacenter.NetworkVm;
import org.cloudbus.cloudsim.network.datacenter.MyParticle;

public class PSO {

	public PSO(double[][] td, double[] et, Datacenter dc)
	{
		workFlowDataTrans = td;
		workFlowTaskExecution = et;
		datacenter = dc;
		setVmData();
	}
	
	public ArrayList<Integer> runPSO()
	{
		System.out.println("workFlowTaskExecution.length: " + workFlowTaskExecution.length);
		// Create a swarm (using 'MyParticle' as sample particle and 'MyFitnessFunction' as fitness function)
		
//		GSA-GA
		Swarm swarm = new SwarmGSA_GA(Swarm.DEFAULT_NUMBER_OF_PARTICLES,
	              new MyParticle(workFlowTaskExecution.length), 
	              new MyFitnessFunction(workFlowDataTrans, workFlowTaskExecution, vmData, vmTransferCost));

		//GSA
//		Swarm swarm = new SwarmGSA(Swarm.DEFAULT_NUMBER_OF_PARTICLES,
//	              new MyParticle(workFlowTaskExecution.length), 
//	              new MyFitnessFunction(workFlowDataTrans, workFlowTaskExecution, vmData, vmTransferCost));
		
//		Myswarm2 swarm = new Myswarm2(Swarm.DEFAULT_NUMBER_OF_PARTICLES,
//	              new MyParticle(workFlowTaskExecution.length), 
//	              new MyFitnessFunction1(workFlowDataTrans, workFlowTaskExecution, vmData, vmTransferCost), 
//	              new MyFitnessFunction2(workFlowDataTrans, workFlowTaskExecution, vmData, vmTransferCost));

		//PSO
//		Swarm swarm = new Swarm(Swarm.DEFAULT_NUMBER_OF_PARTICLES, 
//	               new MyParticle(workFlowTaskExecution.length), 
//	               new MyFitnessFunction(workFlowDataTrans, workFlowTaskExecution, vmData, vmTransferCost));
		
		//System.out.println("workFlowTaskExecution.length: " + workFlowTaskExecution.length);
		// Use neighborhood
		Neighborhood neigh = new Neighborhood1D((Swarm.DEFAULT_NUMBER_OF_PARTICLES-1)/2, true);
		//Neighborhood neigh = new Neighborhood1D(1, true);
		swarm.setNeighborhood(neigh);
		swarm.setNeighborhoodIncrement(0.9);
		
		// Tune swarm's update parameters (if needed)
		swarm.setInertia(0.95);
		swarm.setParticleIncrement(0.8);
		swarm.setGlobalIncrement(0.8);
		
		// Set position (and velocity) constraints. I.e.: where to look for solutions
		swarm.setInertia(0.95);

		System.out.println("setMaxPosition " + vmData.length);
		swarm.setMaxPosition(vmData.length - 1);  //change
		swarm.setMinPosition(0);
		swarm.setMaxMinVelocity(3.5);
		
		int numberOfIterations = NetDatacenterBroker.numberOfIterations;
	
		boolean showGraphics = false;
		
		if (showGraphics) {
			int displayEvery = numberOfIterations / 100 + 1;
			SwarmShow2D ss2d = new SwarmShow2D(swarm, numberOfIterations, displayEvery, true);
			ss2d.run();
		} else {
			// Optimize (and time it)
			for (int i = 0; i < numberOfIterations; i++)
			{
				//STEP1 Initialise if not yet initialised
				//STEP2-3
				swarm.evolve();
				System.out.println(swarm.getBestFitness());
			}
		}
		
		//STEP4 onwards
		// transfer double to int
		double[] bestPosition = swarm.getBestPosition();
		ArrayList<Integer> intBestPosition = new ArrayList<Integer>();
		for(int i = 0; i < bestPosition.length; i++ )
		{
			intBestPosition.add((int)bestPosition[i]);
//			intBestPosition.add(7);
		}
		return intBestPosition;
	}
	
	
//	TODO
	/**
	 * Get VM information(MIPS, execution cost, bandwidth cost) into an array
	 */
	private void setVmData()
	{
	    List<NetworkVm> vmList = datacenter.getVmList();
	    int sizeOfVm = vmList.size();
	    System.out.println("sizeOfVm: " + sizeOfVm);
	    Iterator<NetworkVm> it = vmList.iterator();
	    int count = 0;
	    vmData = new double[sizeOfVm][2];
	    vmTransferCost = new double[sizeOfVm][sizeOfVm];
	    while(it.hasNext())
	    {
	    	NetworkVm tmp = it.next();
	    	vmData[count][0] = tmp.getMips();
	    	vmData[count][1] = tmp.getExecutionCost();
	    	ArrayList<Double> vmBandwidthCost = new ArrayList<Double>();
	    	vmBandwidthCost = tmp.getBandwidthCost();
	    	for (int i = 0; i < vmBandwidthCost.size(); i++)
	    	{
	    		vmTransferCost[count][i] = vmBandwidthCost.get(i);
	    	}
	    	count++;
	    }
	}
	
	
	
	private double[][] vmData;
	private double[][] vmTransferCost;
	private Datacenter datacenter;
	private double[][] workFlowDataTrans;
	private double[] workFlowTaskExecution;
}
