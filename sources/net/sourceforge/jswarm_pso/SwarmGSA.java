package net.sourceforge.jswarm_pso;

import java.awt.Color;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.Iterator;

import org.cloudbus.cloudsim.network.datacenter.MyParticle;

import net.sourceforge.jswarm_pso.FitnessFunction;
import net.sourceforge.jswarm_pso.Neighborhood;
import net.sourceforge.jswarm_pso.Particle;
import net.sourceforge.jswarm_pso.ParticleUpdate;
import net.sourceforge.jswarm_pso.ParticleUpdateSimple;
import net.sourceforge.jswarm_pso.VariablesUpdate;

public class SwarmGSA extends Swarm{

	public static double DEFAULT_GLOBAL_INCREMENT = 0.9;
	public static double DEFAULT_INERTIA = 0.95;
	public static int DEFAULT_NUMBER_OF_PARTICLES = 25;
	public static double DEFAULT_PARTICLE_INCREMENT = 0.9;
	public static double VELOCITY_GRAPH_FACTOR = 10.0;

	/** Best fitness so far (global best) */
	double bestFitness;
	/** Index of best particle so far */
	int bestParticleIndex;
	/** Best position so far (global best) */
	double bestPosition[];
	/** Fitness function for this swarm */
	FitnessFunction fitnessFunction;
	/** Global increment (for velocity update), usually called 'c2' constant */
	double globalIncrement;
	/** Inertia (for velocity update), usually called 'w' constant */
	double inertia;
	/** Maximum position (for each dimension) */
	double maxPosition[];
	/** Maximum Velocity (for each dimension) */
	double maxVelocity[];
	/** Minimum position (for each dimension) */
	double minPosition[];
	/** Minimum Velocity for each dimension. WARNING: Velocity is no in Abs value (so setting minVelocity to 0 is NOT correct!) */
	double minVelocity[];
	/** How many times 'particle.evaluate()' has been called? */
	int numberOfEvaliations;
	/** Number of particles in this swarm */
	int numberOfParticles;
	/** Particle's increment (for velocity update), usually called 'c1' constant */
	double particleIncrement;
	/** Particles in this swarm */
	Particle particles[];
	/** Particle update strategy */
	ParticleUpdate particleUpdate;
	/** A sample particles: Build other particles based on this one */
	Particle sampleParticle;
	/** Variables update */
	VariablesUpdate variablesUpdate;
	/** Neighborhood */
	@SuppressWarnings("unchecked")
	Neighborhood neighborhood;
	/** Neighborhood increment (for velocity update), usually called 'c3' constant */
	double neighborhoodIncrement;
	/** A collection used for 'Iterable' interface */
	ArrayList<Particle> particlesList;
	
	static double[] acceleration;
	static double[][] accelerationArray;
	
	
	double[] allFitness;
	double averageFitness;
	
	//TODO STEP 1
    static int it = 0;
    static int mx_it = 10;
    static double Pm = 0.04;
    static int FEs = 3000;
    
    
	//-------------------------------------------------------------------------
	// Constructors
	//-------------------------------------------------------------------------

	/**
	 * Create a Swarm and set default values
	 * @param numberOfParticles : Number of particles in this swarm (should be greater than 0). 
	 * If unsure about this parameter, try Swarm.DEFAULT_NUMBER_OF_PARTICLES or greater
	 * @param sampleParticle : A particle that is a sample to build all other particles
	 * @param fitnessFunction : Fitness function used to evaluate each particle
	 */
	public SwarmGSA(int numberOfParticles, Particle sampleParticle, FitnessFunction fitnessFunction) {
		super(numberOfParticles, sampleParticle, fitnessFunction);
		super.DEFAULT_NUMBER_OF_PARTICLES=10;
		if (sampleParticle == null) throw new RuntimeException("Sample particle can't be null!");
		if (numberOfParticles <= 0) throw new RuntimeException("Number of particles should be greater than zero.");

		globalIncrement = DEFAULT_GLOBAL_INCREMENT;
		inertia = DEFAULT_INERTIA;
		particleIncrement = DEFAULT_PARTICLE_INCREMENT;
		numberOfEvaliations = 0;
		this.numberOfParticles = numberOfParticles;
		this.sampleParticle = sampleParticle;
		this.fitnessFunction = fitnessFunction;
		bestFitness = Double.NaN;
		bestParticleIndex = -1;

		// Set up particle update strategy (default: ParticleUpdateSimple) 
		particleUpdate = new ParticleUpdateGSA(sampleParticle);

		// Set up variablesUpdate strategy (default: VariablesUpdate)
		variablesUpdate = new VariablesUpdate();

		neighborhood = null;
		neighborhoodIncrement = 0.0;
		particlesList = null;
	}

	//-------------------------------------------------------------------------
	// Methods
	//-------------------------------------------------------------------------

	public double getBestFitness() {
		return bestFitness;
	}

	public Particle getBestParticle() {
		return particles[bestParticleIndex];
	}

	public int getBestParticleIndex() {
		return bestParticleIndex;
	}

	public double[] getBestPosition() {
		return bestPosition;
	}

	public FitnessFunction getFitnessFunction() {
		return fitnessFunction;
	}

	public double getGlobalIncrement() {
		return globalIncrement;
	}

	public double getInertia() {
		return inertia;
	}

	public double[] getMaxPosition() {
		return maxPosition;
	}

	public double[] getMaxVelocity() {
		return maxVelocity;
	}

	public double[] getMinPosition() {
		return minPosition;
	}

	public double[] getMinVelocity() {
		return minVelocity;
	}

	@SuppressWarnings("unchecked")
	public Neighborhood getNeighborhood() {
		return neighborhood;
	}

	/**
	 * Return the best position in the neighborhood
	 * Note: If neighborhood is not defined (i.e. neighborhood is null) then 'particle' is returned 
	 * so that it doesn't influence in particle update.
	 * 
	 * @param particle
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public double[] getNeighborhoodBestPosition(Particle particle) {
		if (neighborhood == null) return particle.getPosition();
		double d[] = neighborhood.getBestPosition(particle);
		if (d == null) return particle.getPosition();
		return d;
	}

	public double getNeighborhoodIncrement() {
		return neighborhoodIncrement;
	}

	public int getNumberOfEvaliations() {
		return numberOfEvaliations;
	}

	public int getNumberOfParticles() {
		return numberOfParticles;
	}

	public Particle getParticle(int i) {
		return particles[i];
	}

	public double getParticleIncrement() {
		return particleIncrement;
	}

	public Particle[] getParticles() {
		return particles;
	}

	public ParticleUpdate getParticleUpdate() {
		return particleUpdate;
	}

	public Particle getSampleParticle() {
		return sampleParticle;
	}

	public VariablesUpdate getVariablesUpdate() {
		return variablesUpdate;
	}

	public static double[] getAcceleration() {
		return acceleration;
	}
	
	/**
	 * Make an iteration: 
	 * 	- evaluates the swarm 
	 * 	- updates positions and velocities
	 * 	- applies positions and velocities constraints 
	 */
	public void evolve() {
		// Initialize (if not already done)
		if (particles == null) init();		//STEP1

		evaluate(); 						// Evaluate particles STEP2
//		GSA_GA();							// STEP3
		
		update2(); // Update positions and velocities

		variablesUpdate.update(this);
	}
	
	/**
	 * Initialize every particle
	 * Warning: maxPosition[], minPosition[], maxVelocity[], minVelocity[] must be initialized and setted
	 * TODO STEP1
	 */
	public void init() {
		// Init particles
		particles = new Particle[numberOfParticles];

		// Check constraints (they will be used to initialize particles)
		if (maxPosition == null) throw new RuntimeException("maxPosition array is null!");
		if (minPosition == null) throw new RuntimeException("maxPosition array is null!");
		if (maxVelocity == null) {
			// Default maxVelocity[]
			int dim = sampleParticle.getDimension();
			maxVelocity = new double[dim];
			for (int i = 0; i < dim; i++)
				maxVelocity[i] = (maxPosition[i] - minPosition[i]) / 2.0;
		}
		if (minVelocity == null) {
			// Default minVelocity[]
			int dim = sampleParticle.getDimension();
			minVelocity = new double[dim];
			for (int i = 0; i < dim; i++)
				minVelocity[i] = -maxVelocity[i];
		}

		// Init each particle RANDOMLY
		System.out.println("numberOfParticles: " + numberOfParticles);
		for (int i = 0; i < numberOfParticles; i++) {
			System.out.println("maxPosition: " + maxPosition[0]);
			System.out.println("minPosition: " + minPosition[0]);
			particles[i] = (Particle) sampleParticle.selfFactory(); // Create a new particles (using 'sampleParticle' as reference)
			particles[i].init(maxPosition, minPosition, maxVelocity, minVelocity); // Initialize it
		}

		// Init neighborhood
		if (neighborhood != null) neighborhood.init(this);
	}

	/**
	 * Evaluate fitness function for every particle 
	 * Warning: particles[] must be initialized and fitnessFunction must be set
	 * TODO STEP2
	 */
	public void evaluate() {
		if (particles == null) throw new RuntimeException("No particles in this swarm! May be you need to call Swarm.init() method");
		if (fitnessFunction == null) throw new RuntimeException("No fitness function in this swarm! May be you need to call Swarm.setFitnessFunction() method");

		// Initialize
		if (Double.isNaN(bestFitness)) {
			bestFitness = (fitnessFunction.isMaximize() ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY);
			bestParticleIndex = -1;
		}

		System.out.println(fitnessFunction.isMaximize());
		
		allFitness = new double[particles.length];
		averageFitness = 0;
		//---
		// Evaluate each particle (and find the 'best' one)
		//---
		for (int i = 0; i < particles.length; i++) {
			// Evaluate particle
			System.out.println("particle" + particles[i].getPosition()[0]);
			System.out.println("particles.length" + particles.length);
			double fit = fitnessFunction.evaluate(particles[i]);
			allFitness[i] = fit;
			averageFitness += fit;
			
			numberOfEvaliations++; // Update counter

			// Update 'best global' position
			if (fitnessFunction.isBetterThan(bestFitness, fit)) {
				bestFitness = fit; // Copy best fitness, index, and position vector
				bestParticleIndex = i;
				if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
				particles[bestParticleIndex].copyPosition(bestPosition);
			}

			// Update 'best neighborhood' 
			if (neighborhood != null) {
				neighborhood.update(this, particles[i]);
			}
		}
		System.out.println("BEST FITNESS:  " + bestFitness + "\nBEST POSITION:  ");
		for (int i=0; i<bestPosition.length; i++ ){ 
			System.out.print((int)bestPosition[i] + " ");
		}
		System.out.println(bestParticleIndex+ "");
		averageFitness /= particles.length;
	}
	
	/**
	 * TODO STEP3
	 */
	public void GSA_GA() {
		
		//create random number
		double rand = Math.random();
		
		if(rand>standard_deviation() || rand<it/mx_it) {
			// call GA
			GAalgo();		// replace by GAalgo() STEP3a
			System.out.println("GA called");
		}
		else {
			// call GSA
			update2();		// STEP3b
			System.out.println("GSA called");
		}
		
	}
	
	public double standard_deviation() {
		double sd = 0;
		for (int i=0; i<allFitness.length;i++)
		{
		    sd = sd + Math.pow(allFitness[i] - averageFitness, 2);
		}
		return sd;
	}

	/**
	 * Genetic Algorithm
	 * TODO STEP3a
	 */
	public void GAalgo() {
		
		// Roulette Selection of parents
		double totalFitness = averageFitness * particles.length;
		double[] prob = new double[particles.length+1];
		for(int i=1; i<particles.length; i++)
			prob[i] = (allFitness[i]/totalFitness);
		
		ArrayList<Integer> toMutate = new ArrayList<>();
		for(int i=0; i<particles.length; i++) {
			if(Pm>prob[i])
				toMutate.add(i);
		}
		
		double k0 = 1.5;		// mutation constant
		
		// Mutation TODO
		for(int i=0; i<toMutate.size(); i++) {
			int pos = toMutate.get(i);
			allFitness[pos] *= k0;
			
			// Update best position based on GA TODO
			double fit = allFitness[pos];
			if (fitnessFunction.isBetterThan(bestFitness, fit)) {
				bestFitness = fit; // Copy best fitness, index, and position vector
				bestParticleIndex = i;
				if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
				particles[bestParticleIndex].copyPosition(bestPosition);
			}
		}
		
		update2();
		
	}
	
	/**
	 * Iterate over all particles
	 */
	public Iterator<Particle> iterator() {
		if (particlesList == null) {
			particlesList = new ArrayList<Particle>(particles.length);
			for (int i = 0; i < particles.length; i++)
				particlesList.add(particles[i]);
		}

		return particlesList.iterator();
	}

	public void setBestParticleIndex(int bestParticle) {
		bestParticleIndex = bestParticle;
	}

	public void setBestPosition(double[] bestPosition) {
		this.bestPosition = bestPosition;
	}

	public void setFitnessFunction(FitnessFunction fitnessFunction) {
		this.fitnessFunction = fitnessFunction;
	}

	public void setGlobalIncrement(double globalIncrement) {
		this.globalIncrement = globalIncrement;
	}

	public void setInertia(double inertia) {
		this.inertia = inertia;
	}

	/**
	 * Sets every maxVelocity[] and minVelocity[] to 'maxVelocity' and '-maxVelocity' respectively
	 * @param maxVelocity
	 */
	public void setMaxMinVelocity(double maxVelocity) {
		if (sampleParticle == null) throw new RuntimeException("Need to set sample particle before calling this method (use Swarm.setSampleParticle() method)");
		int dim = sampleParticle.getDimension();
		this.maxVelocity = new double[dim];
		minVelocity = new double[dim];
		for (int i = 0; i < dim; i++) {
			this.maxVelocity[i] = maxVelocity;
			minVelocity[i] = -maxVelocity;
		}
	}

	/**
	 * Sets every maxPosition[] to 'maxPosition'
	 * @param maxPosition
	 */
	public void setMaxPosition(double maxPosition) {
		if (sampleParticle == null) throw new RuntimeException("Need to set sample particle before calling this method (use Swarm.setSampleParticle() method)");
		int dim = sampleParticle.getDimension();
		this.maxPosition = new double[dim];
		for (int i = 0; i < dim; i++)
			this.maxPosition[i] = maxPosition;
	}

	public void setMaxPosition(double[] maxPosition) {
		this.maxPosition = maxPosition;
	}

	public void setMaxVelocity(double[] maxVelocity) {
		this.maxVelocity = maxVelocity;
	}

	/**
	 * Sets every minPosition[] to 'minPosition'
	 * @param minPosition
	 */
	public void setMinPosition(double minPosition) {
		if (sampleParticle == null) throw new RuntimeException("Need to set sample particle before calling this method (use Swarm.setSampleParticle() method)");
		int dim = sampleParticle.getDimension();
		this.minPosition = new double[dim];
		for (int i = 0; i < dim; i++)
			this.minPosition[i] = minPosition;
	}

	public void setMinPosition(double[] minPosition) {
		this.minPosition = minPosition;
	}

	public void setMinVelocity(double minVelocity[]) {
		this.minVelocity = minVelocity;
	}

//	@SuppressWarnings("unchecked")
	public void setNeighborhood(Neighborhood neighborhood) {
		this.neighborhood = neighborhood;
	}

	public void setNeighborhoodIncrement(double neighborhoodIncrement) {
		this.neighborhoodIncrement = neighborhoodIncrement;
	}

	public void setNumberOfEvaliations(int numberOfEvaliations) {
		this.numberOfEvaliations = numberOfEvaliations;
	}

	public void setNumberOfParticles(int numberOfParticles) {
		this.numberOfParticles = numberOfParticles;
	}

	public void setParticleIncrement(double particleIncrement) {
		this.particleIncrement = particleIncrement;
	}

	public void setParticles(Particle[] particle) {
		particles = particle;
		particlesList = null;
	}

	public void setParticleUpdate(ParticleUpdate particleUpdate) {
		this.particleUpdate = (ParticleUpdateGSA) particleUpdate;
	}

	public void setSampleParticle(Particle sampleParticle) {
		this.sampleParticle = sampleParticle;
	}

	public void setVariablesUpdate(VariablesUpdate variablesUpdate) {
		this.variablesUpdate = variablesUpdate;
	}

	/**
	 * Show a swarm in a graph 
	 * @param graphics : Grapics object
	 * @param foreground : foreground color
	 * @param width : graphic's width
	 * @param height : graphic's height
	 * @param dim0 : Dimention to show ('x' axis)
	 * @param dim1 : Dimention to show ('y' axis)
	 * @param showVelocity : Show velocity tails?
	 */
	public void show(Graphics graphics, Color foreground, int width, int height, int dim0, int dim1, boolean showVelocity) {
		graphics.setColor(foreground);

		if (particles != null) {
			double scalePosW = width / (maxPosition[dim0] - minPosition[dim0]);
			double scalePosH = height / (maxPosition[dim1] - minPosition[dim1]);
			double minPosW = minPosition[dim0];
			double minPosH = minPosition[dim1];

			double scaleVelW = width / (VELOCITY_GRAPH_FACTOR * (maxVelocity[dim0] - minVelocity[dim0]));
			double scaleVelH = height / (VELOCITY_GRAPH_FACTOR * (maxVelocity[dim1] - minVelocity[dim1]));
			double minVelW = minVelocity[dim0] + (maxVelocity[dim0] - minVelocity[dim0]) / 2;
			double minVelH = minVelocity[dim1] + (maxVelocity[dim1] - minVelocity[dim1]) / 2;

			for (int i = 0; i < particles.length; i++) {
				int vx, vy, x, y;
				double pos[] = particles[i].getPosition();
				double vel[] = particles[i].getVelocity();
				x = (int) (scalePosW * (pos[dim0] - minPosW));
				y = height - (int) (scalePosH * (pos[dim1] - minPosH));
				graphics.drawRect(x - 1, y - 1, 3, 3);
				if (showVelocity) {
					vx = (int) (scaleVelW * (vel[dim0] - minVelW));
					vy = (int) (scaleVelH * (vel[dim1] - minVelH));
					graphics.drawLine(x, y, x + vx, y + vy);
				}
			}
		}
	}

	/** Swarm size (number of particles) */
	public int size() {
		return particles.length;
	}

	/** Printable string */
	@Override
	public String toString() {
		String str = "";

		if (particles != null) str += "Swarm size: " + particles.length + "\n";

		if ((minPosition != null) && (maxPosition != null)) {
			str += "Position ranges:\t";
			for (int i = 0; i < maxPosition.length; i++)
				str += "[" + minPosition[i] + ", " + maxPosition[i] + "]\t";
		}

		if ((minVelocity != null) && (maxVelocity != null)) {
			str += "\nVelocity ranges:\t";
			for (int i = 0; i < maxVelocity.length; i++)
				str += "[" + minVelocity[i] + ", " + maxVelocity[i] + "]\t";
		}

		if (sampleParticle != null) str += "\nSample particle: " + sampleParticle;

		if (particles != null) {
			str += "\nParticles:";
			for (int i = 0; i < particles.length; i++) {
				str += "\n\tParticle: " + i + "\t";
				str += particles[i].toString();
			}
		}
		str += "\n";

		return str;
	}

	/**
	 * Return a string with some (very basic) statistics 
	 * @return A string
	 */
	public String toStringStats() {
		String stats = "";
		if (!Double.isNaN(bestFitness)) {
			stats += "Best fitness: " + bestFitness + "\nBest position: \t[";
			for (int i = 0; i < bestPosition.length; i++)
				stats += bestPosition[i] + (i < (bestPosition.length - 1) ? ", " : "");
			stats += "]\nNumber of evaluations: " + numberOfEvaliations + "\n";
		}
		return stats;
	}

	/**
	 * Update every particle's position and velocity, also apply position and velocity constraints (if any)
	 * Warning: Particles must be already evaluated
	 */
	public void update() {
		it++;
		
		particleUpdate.begin(this);
		
		int dimension=particles.length;
		double worstFitness=1000000000,bestFitness=-1000000000;
		for(int i=0;i < dimension; i++){
			bestFitness=Math.min(bestFitness, particles[i].getFitness());
			worstFitness=Math.max(worstFitness, particles[i].getFitness());
		}
		double imass[] = new double[dimension];
		double gmass[] = new double[dimension];
		double gravity[] = new double[dimension];
		double force[] = new double[dimension];
		double tin=0.0;
		double acc[] = new double[dimension];

//		for(int i=0;i < dimension;i++)
//			for(int j=0; j<particles[0].position.length;j++){
//				particles[i].position[j]=particles[i].velocity[j]=0;
//				//System.out.println(i +" " + j);
//			}
		double ttin = 0.0;
		for(int i=0 ;i < dimension ;i++){
			imass[i]=(particles[i].getFitness()-worstFitness)/(bestFitness-worstFitness);
			ttin+=imass[i]; 
		}
		tin = ttin;
	   for(int i=0; i< dimension; i++){
		   gmass[i]=imass[i]/tin;
	   }
	  // System.out.println("dimension " + dimension);
	   for(int i = 0; i < dimension; i++){
		   force[i] = 0.0;
		   double[] pos = particles[i].getPosition();
		   for(int j = 0; j < pos.length; j++){
			   if(i == j)
				   continue;
			   //System.out.println(i+" "+j);
				   //System.out.println(dimension);
			 //  System.out.println(i+" "+j+" "+particles[i].position[j]);
			   double temp = Math.random() * ((getG(it,mx_it) * gmass[i] * gmass[j])*(pos[j])) / (euclidDis(i,j));			//me removed *Math.sqrt((gmass[i]+gmass[j])/2)/getG(it,mx_it) from denominator
			   force[i] += temp;
		   }
		   acc[i] = force[i] / gmass[i];		 //me had to uncomment this
	   }
	   
	   for(int i = 0; i < dimension; i++){
		   double[] vel = particles[i].getVelocity();
		   double[] pos = particles[i].getPosition();
		   for(int j=0;j < pos.length; j++){
			   vel[j] = Math.random() * vel[j] + acc[i];
			   pos[j] += vel[j];
		  // System.out.println(i+" "+j+" "+particles[i].velocity[j]);
		   }
		   
		   //System.out.println(particles[i].getDimension());
		   System.out.println("maxPosition: " + maxPosition[0]);
		   particles[i].applyConstraints(minPosition, maxPosition, minVelocity, maxVelocity);
	   }
	  // System.out.println("dimension " + dimension);
	   particleUpdate.end(this);
	}
	
	public void update2() {

		it++;
		
		// Initialize a particle update iteration
		particleUpdate.begin(this);

		setAccelerationArray();
		
		// For each particle...
		for (int i = 0; i < particles.length; i++) {
			// Update particle's position and speed
//			acceleration = accelerationArray[i];
			acceleration = accelerationArray[i];
			particleUpdate.update(this, particles[i]);

			// Apply position and velocity constraints
			particles[i].applyConstraints(minPosition, maxPosition, minVelocity, maxVelocity);
		}

		// Finish a particle update iteration
		particleUpdate.end(this);
	}

	public void setAccelerationArray() {

		int N=particles.length;
		int dimension = particles[0].getPosition().length;
		
		double worstFitness=1000000000,bestFitness=-1000000000;
		for(int i=0;i < N; i++){
			System.out.println(particles[i].getFitness() + ", " + allFitness[i]);
//			bestFitness=Math.max(bestFitness, particles[i].getFitness());
//			worstFitness=Math.min(worstFitness, particles[i].getFitness());
			bestFitness=Math.max(bestFitness, allFitness[i]);
			worstFitness=Math.min(worstFitness, allFitness[i]);
		}
		System.out.println("");
		System.out.println(bestFitness);
		System.out.println(worstFitness);
		double imass[] = new double[N];
		double gmass[] = new double[N];
		double gravity[] = new double[N];
		double force[][] = new double[N][dimension];
		double tin=0.0;
		double acc[][] = new double[N][dimension];
		

//		for(int i=0;i < dimension;i++)
//			for(int j=0; j<particles[0].position.length;j++){
//				particles[i].position[j]=particles[i].velocity[j]=0;
//			}
		double ttin = 0.0;
		for(int i=0 ;i < N ;i++){
//			imass[i]=(particles[i].getFitness())/(bestFitness);
			imass[i]=(allFitness[i] - worstFitness)/(bestFitness - worstFitness);
			ttin+=imass[i];
		}
		
		for(int i=0; i<N; i++){
			imass[i] = imass[i]/(ttin * bestFitness);
			tin += imass[i];
			
			System.out.print(imass[i] + ", ");
		}
		System.out.println("");
	   for(int i=0; i< N; i++){
		   gmass[i]=imass[i]/tin;
			System.out.print(gmass[i] + ", ");
	   }
		System.out.println("");
	   double total_acc = 0.0;
	  // System.out.println("dimension " + dimension);
	   for(int i = 0; i < N; i++){
//		   force[i] = 0.0;
		   double[] pos = particles[i].getPosition();
		   double[] vel = particles[i].getVelocity(); 
		   for(int k=0 ; k< dimension; k++) {
			   force[i][k] = 0;
			   acc[i][k] = 0;
			   for(int j = 0; j < N; j++){
				   
				   double[] pos2 = particles[j].getPosition();
				   double[] vel2 = particles[j].getVelocity();
				   
				   double euclid = euclidDis(i,j);

//				   System.out.println(euclid);
//				   System.out.println(i+" "+j);
					   //System.out.println(dimension);
				 //  System.out.println(i+" "+j+" "+particles[i].position[j]);
				   double temp = Math.random() * ( getG(it,mx_it) *gmass[i]*gmass[j]* direction(pos2[k], pos[k]) / euclid);
				   force[i][k] += temp;
				   
				   acc[i][k] += force[i][k];
//				   acc[i][k] += (pos2[k] - pos[k]);
			   }
//			   System.out.print(acc[i][k] + ",  ");
		   }
//		   System.out.println("");
	   }
//	   for(int j=0; j<dimension; j++) {
//		   double tot_acc = 0.0;
//		   for(int i=0; i<N; i++) {
//			   tot_acc += modulus(acc[i][j]);
//		   }
//		   for(int i=0; i<N; i++) {
//			   acc[i][j] /= tot_acc;
//			   System.out.print(acc[i][j] + ",  ");
//		   }
//		   System.out.println("");
//	   }

	   SwarmGSA.accelerationArray = acc;
	}
	
	double euclidDis(int u,int v){
		double eps = 1;
		double dis=0.0;
		double[] pos1 = particles[u].getPosition();
		double[] pos2 = particles[v].getPosition();
		for(int i=0;i<pos1.length;i++){
			dis=dis+(pos1[i]-pos2[i])*(pos1[i]-pos2[i]);
		}
//		dis=Math.sqrt(dis);
		if(dis!=0)
			return dis;
		else
			return eps;
	}
	
	
	public double getG(int it,int mx_it){
		double alfa=20,g0=10;
		double G=g0*(Math.exp(-1*alfa*it/mx_it));
//		double G=g0*(1-it/(alfa*mx_it));
		return G;
	}

	public double modulus(double x) {
		if(x >= 0)
			return x;
		else
			return -1 * x;
	}
	public double direction(double x, double y) {
		if(x != y)
			return (x-y)/modulus(x-y);
		else
			return 1;
	}
	
}
//450053.85111111094
//36388.33
//TODO cloudlet parameters
//cloudletLength
//cloudletFileSize
//numberOfPes
