package org.cloudbus.cloudsim.network.datacenter;

import net.sourceforge.jswarm_pso.FitnessFunction;

public class MyFitnessFunction2 extends FitnessFunction{
	public MyFitnessFunction2(double[][] td, double[] et, double[][] vd, double[][] vt){
		super(false);
		workFlowDataTrans = td;
		workFlowTaskExcution = et;
		vmData = vd;
		vmTransferCost = vt;
	}
	
	public double evaluate(double position[]) {
		double fitnessValue = 0;
		int[] intPosition = new int[position.length];
		//System.out.println("position.length: " + position.length);
		for(int i = 0; i < position.length; i++ )
		{
			intPosition[i] = (int)position[i];
		}
		
		double[] vmCost = new double[vmData.length];
		double [] vmtime =new double[vmData.length];
		// calculate each task's total cost, and add to the vm's cost the task assigned to
//		System.out.println("abcd");
		for(int i = 0; i < position.length; i++ )
		{
			int vmNum = intPosition[i];
			vmtime[vmNum]+=workFlowTaskExcution[i] / vmData[vmNum][0];
			System.out.print(intPosition[i]+"+");

			double taskCost = workFlowTaskExcution[i] / vmData[vmNum][0] * vmData[vmNum][1]; //add execution cost
			//System.out.println("workFlowTaskExcution[i]: " + workFlowTaskExcution[i]);
			//System.out.println("vmData[vmNum][0]" + vmData[vmNum][0]);
			//System.out.println("vmData[vmNum][1]" + vmData[vmNum][1]);
			double taskDataTransferCost = 0;
			for(int j = 0; j < workFlowDataTrans[i].length && j < intPosition.length; j++)
			{
//				System.out.println(i+", "+j);
				if(workFlowDataTrans[i][j] != 0)
				{
//					System.out.println(workFlowDataTrans[i][j]);
//					System.out.println(vmTransferCost[vmNum][intPosition[j]]);
					taskDataTransferCost += workFlowDataTrans[i][j] * vmTransferCost[vmNum][intPosition[j]]; //add the task output file transfer cost
				}
				
			}
			//System.out.println("taskCost: " + taskCost);
			//System.out.println("taskDataTransferCost: " + taskDataTransferCost);
			vmCost[vmNum] += taskDataTransferCost + taskCost;
			//System.out.println("vmCost[vmNum]: " + vmCost[vmNum]);
		}
		double alpha=1.0;
		
		//System.out.println("fitnessValue: " + nfit);
		
		fitnessValue = max(vmCost);
		System.out.println();
		System.out.println("fitnessValue1: " + fitnessValue);
		System.out.println("fitnessValue2: " + max(vmtime));
		System.out.println("fitnessValue3: " + alpha*fitnessValue+(1-alpha)*(max(vmtime)));
		//System.out.println("fitnessValue: " +alpha*fitnessValue+(1-alpha)*(max(vmtime)));
		
		return max(vmtime);
		
		//return fitnessValue;
	}
	
	public double[] evaluate1(double position[]) {
		double returnValue[] = new double[3];
		double fitnessValue = 0;
		int[] intPosition = new int[position.length];
		//System.out.println("position.length: " + position.length);
		for(int i = 0; i < position.length; i++ )
		{
			intPosition[i] = (int)position[i];
		}
		
		double[] vmCost = new double[vmData.length];
		double [] vmtime =new double[vmData.length];
		// calculate each task's total cost, and add to the vm's cost the task assigned to
		for(int i = 0; i < position.length; i++ )
		{
			int vmNum = intPosition[i];
			vmtime[vmNum]+=workFlowTaskExcution[i] / vmData[vmNum][0];
			System.out.print(intPosition[i]+"+");

			double taskCost = workFlowTaskExcution[i] / vmData[vmNum][0] * vmData[vmNum][1]; //add execution cost
			//System.out.println("workFlowTaskExcution[i]: " + workFlowTaskExcution[i]);
			//System.out.println("vmData[vmNum][0]" + vmData[vmNum][0]);
			//System.out.println("vmData[vmNum][1]" + vmData[vmNum][1]);
			double taskDataTransferCost = 0;
			for(int j = 0; j < workFlowDataTrans[i].length; j++)
			{
				if(workFlowDataTrans[i][j] != 0)
				{
					taskDataTransferCost += workFlowDataTrans[i][j] * vmTransferCost[vmNum][intPosition[j]]; //add the task output file transfer cost
				}
				
			}
			//System.out.println("taskCost: " + taskCost);
			//System.out.println("taskDataTransferCost: " + taskDataTransferCost);
			vmCost[vmNum] += taskDataTransferCost + taskCost;
			//System.out.println("vmCost[vmNum]: " + vmCost[vmNum]);
		}
		double alpha=1.0;
		
		//System.out.println("fitnessValue: " + nfit);
		
		fitnessValue = max(vmCost);
//		System.out.println("fitnessValue1: " + fitnessValue);
//		System.out.println("fitnessValue2: " + max(vmtime));
		//System.out.println("fitnessValue: " +alpha*fitnessValue+(1-alpha)*(max(vmtime)));
		returnValue[0] = alpha*fitnessValue+(1-alpha)*(max(vmtime));
		returnValue[1] = fitnessValue;
		returnValue[2] = max(vmtime);
		
		return returnValue;
		
		//return fitnessValue;
	}
	
	public double getfitness1(double position[]) {
		
		
		return 0;
	}
	
	/**
	 * get max value in a array
	 * @param list input array
	 * @return the max vlaue
	 */
	private double max(double[] list)
	{
		double max = 0;
		for(int i = 0; i < list.length; i++ )
		{
			if (list[i] > max)
			{
				max = list[i];
			}	
		}
		return max;
	}
	
	private double[][] workFlowDataTrans;
	private double[] workFlowTaskExcution;
	private double[][] vmData;
	private double[][] vmTransferCost;
}
