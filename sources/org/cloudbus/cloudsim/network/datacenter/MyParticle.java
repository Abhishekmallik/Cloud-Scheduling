package org.cloudbus.cloudsim.network.datacenter;

import net.sourceforge.jswarm_pso.Particle;

public class MyParticle extends Particle
{
	public MyParticle(int d)
	{
		super(d);
		//System.out.println("d:" + d);
	}
	
	public MyParticle()
	{
		super(20);
		//System.out.println("10");
		//need to change the number of tasks here
		//double acc[]=new double[10],force[]=new double[10],amass[]=new double[10],pmass[]=new double[10];
	}
}
