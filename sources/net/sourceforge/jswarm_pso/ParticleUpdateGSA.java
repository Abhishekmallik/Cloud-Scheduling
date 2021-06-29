
package net.sourceforge.jswarm_pso;

/**
 * Particle update strategy
 * 
 * Every Swarm.evolve() itereation the following methods are called
 * 		- begin(Swarm) : Once at the begining of each iteration
 * 		- update(Swarm,Particle) : Once for each particle
 * 		- end(Swarm) : Once at the end of each iteration
 * 
 * @author Pablo Cingolani <pcingola@users.sourceforge.net>
 */
public class ParticleUpdateGSA extends ParticleUpdate {

	/** Random vector for local update */
	double rlocal[];
	/** Random vector for global update */
	double rglobal[];
	/** Random vector for neighborhood update */
	double rneighborhood[];

	/**
	 * Constructor 
	 * @param particle : Sample of particles that will be updated later
	 */
	public ParticleUpdateGSA(Particle particle) {
		super(particle);
		rlocal = new double[particle.getDimension()];
		rglobal = new double[particle.getDimension()];
		rneighborhood = new double[particle.getDimension()];
	}

	/** 
	 * This method is called at the begining of each iteration
	 * Initialize random vectors use for local and global updates (rlocal[] and rother[])
	 */
	@Override
	public void begin(Swarm swarm) {
		int i, dim = swarm.getSampleParticle().getDimension();
		for (i = 0; i < dim; i++) {
			rlocal[i] = Math.random();
			rglobal[i] = Math.random();
			rneighborhood[i] = Math.random();
		}
	}

	/** This method is called at the end of each iteration */
	@Override
	public void end(Swarm swarm) {
	}

	/** Update particle's velocity and position */
	@Override
	public void update(Swarm swarm, Particle particle) {
		
		double[] acc = SwarmGSA.getAcceleration();
//		System.out.println(acc);
		
		double position[] = particle.getPosition();
		double velocity[] = particle.getVelocity();
		double fitness = particle.getFitness();
		
		double globalBestPosition[] = swarm.getBestPosition();
		double particleBestPosition[] = particle.getBestPosition();
		double neighBestPosition[] = swarm.getNeighborhoodBestPosition(particle);

		// Update velocity and position
		for (int i = 0; i < position.length; i++) {
			// Update velocity
//			velocity[i] = swarm.getInertia() * velocity[i] // Inertia
//					+ rlocal[i] * swarm.getParticleIncrement() * (particleBestPosition[i] - position[i]) // Local best
//					+ rneighborhood[i] * swarm.getNeighborhoodIncrement() * (neighBestPosition[i] - position[i]) // Neighborhood best					
//					+ rglobal[i] * swarm.getGlobalIncrement() * (globalBestPosition[i] - position[i]); // Global best
			// Update position
//			System.out.print(swarm.getInertia() * velocity[i] 
//					+ rlocal[i] * swarm.getParticleIncrement() * (particleBestPosition[i] - position[i]) 
//					+ rneighborhood[i] * swarm.getNeighborhoodIncrement() * (neighBestPosition[i] - position[i]) // Neighborhood best					
//					+ rglobal[i] * swarm.getGlobalIncrement() * (globalBestPosition[i] - position[i]) // Global best
//					
//			 + ",  ");
			System.out.print(acc[i] + ", ");
			velocity[i] = velocity[i] + acc[i];
			position[i] += velocity[i];
		}
		System.out.println("");
		
	}
}
