package net.sourceforge.jswarm_pso;

import java.awt.Color;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

import org.cloudbus.cloudsim.network.datacenter.MyParticle;
import org.cloudbus.cloudsim.network.datacenter.NetDatacenterBroker;

import net.sourceforge.jswarm_pso.FitnessFunction;
import net.sourceforge.jswarm_pso.Neighborhood;
import net.sourceforge.jswarm_pso.Particle;
import net.sourceforge.jswarm_pso.ParticleUpdate;
import net.sourceforge.jswarm_pso.ParticleUpdateSimple;
import net.sourceforge.jswarm_pso.VariablesUpdate;

public class SwarmGSA_GA extends Swarm{

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
	public SwarmGSA_GA(int numberOfParticles, Particle sampleParticle, FitnessFunction fitnessFunction) {
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
		particleUpdate = new ParticleUpdateGSA_GA(sampleParticle);

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

		for(int i=0; i<particles.length; i++) {
			for(int j=0; j<particles[i].position.length; j++) {
				System.out.print(particles[i].position[j] + ",");
			}
			System.out.println("");
		}
		for(int i=0; i<particles.length; i++) {
			for(int j=0; j<particles[i].velocity.length; j++) {
				System.out.print(particles[i].velocity[j] + ",");
			}
			System.out.println("");
		}
		
		evaluate(); 						// Evaluate particles STEP2
		GSA_GA();							// STEP3
		
//		update2(); // Update positions and velocities

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

		int di = sampleParticle.getDimension();
		System.out.println("DIMENSION = " + di);
		// Check constraints (they will be used to initialize particles)
		if (maxPosition == null) throw new RuntimeException("maxPosition array is null!");
		if (minPosition == null) throw new RuntimeException("maxPosition array is null!");
		if (maxVelocity == null) {
			// Default maxVelocity[]
			int dim = sampleParticle.getDimension();
			System.out.println("DIMENSION = " + dim);
			maxVelocity = new double[dim];
			for (int i = 0; i < dim; i++) {
				maxVelocity[i] = (maxPosition[i] - minPosition[i]) / 2.0;

			}
		}
		if (minVelocity == null) {
			// Default minVelocity[]
			int dim = sampleParticle.getDimension();
			minVelocity = new double[dim];
			for (int i = 0; i < dim; i++) {
				minVelocity[i] = -maxVelocity[i];
				System.out.println("maxVelocity: " + maxVelocity[i]);
				System.out.println("minVelocity: " + minVelocity[i]);
			}
		}

		// Init each particle RANDOMLY
		System.out.println("numberOfParticles: " + numberOfParticles);
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i] = (Particle) sampleParticle.selfFactory(); // Create a new particles (using 'sampleParticle' as reference)
			System.out.println(particles[i].position.length);
			particles[i].init(maxPosition, minPosition, maxVelocity, minVelocity); // Initialize it
//			particles[i].position = particle_position()[i];
//			particles[i].velocity = particle_velocity()[i];
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
//		double rand = Math.random();
		double rand = 0.5;
		
		if(rand>standard_deviation() || rand>it/mx_it) {
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
		double[] prob = new double[particles.length];
		for(int i=0; i<particles.length; i++)
			prob[i] = (allFitness[i]/totalFitness);

		ArrayList<Integer> toMutate = new ArrayList<>();
		ArrayList<Integer> notToMutate = new ArrayList<>();
		for(int i=0; i<particles.length; i++) {
			if(Pm<prob[i])
				toMutate.add(i);
			else
				notToMutate.add(i);
		}
		
		// Cross-over TODO
		crossover1(toMutate, notToMutate);
		
		// Mutation TODO
		mutate(toMutate, notToMutate);
		
//		for(int i=0; i<toMutate.size(); i++) {
//			int pos = toMutate.get(i);
//			allFitness[pos] *= k0;
//			
//			// Update best position based on GA TODO
//			double fit = allFitness[pos];
//			if (fitnessFunction.isBetterThan(bestFitness, fit)) {
//				bestFitness = fit; // Copy best fitness, index, and position vector
//				bestParticleIndex = i;
//				if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
//				particles[bestParticleIndex].copyPosition(bestPosition);
//			}
//		}
		
		update2();
		
	}

	private void crossover1(ArrayList<Integer> toMutate, ArrayList<Integer> notToMutate) {

		int dimen = particles[0].position.length;
		int mut_size = notToMutate.size();
		
		for(int i=0; i<toMutate.size(); i++) {
			int rand = (int)Math.random()*mut_size;
			
			//parents' position
			double[] pos_parent1 = particles[notToMutate.get(rand)].position;		//better parent
			double[] pos_parent2 = particles[toMutate.get(i)].position;				//worse  parent
			
			//offsprings' position
			double[] pos_child1 = new double[dimen];
			double[] pos_child2 = new double[dimen];
			
			int half = dimen/2;
			for(int j=0; j<half; j++) {
				pos_child1[j] = pos_parent1[j];
				pos_child1[j + half] = pos_parent2[j + half];
				
				pos_child2[j] = pos_parent2[j];
				pos_child2[j + half] = pos_parent1[j + half];
			}
			if(dimen%2 != 0) {
				pos_child1[dimen-1] = pos_parent2[dimen-1]; 
				pos_child2[dimen-1] = pos_parent1[dimen-1];
			}

			double fitness_child1 = fitnessFunction.evaluate(pos_child1);
			double fitness_child2 = fitnessFunction.evaluate(pos_child2);
			double parent_fitness = allFitness[toMutate.get(i)];

			if (fitnessFunction.isBetterThan(parent_fitness, fitness_child1) && fitnessFunction.isBetterThan(fitness_child2, fitness_child1)) {
				particles[toMutate.get(i)].position = pos_child1;

				if (fitnessFunction.isBetterThan(bestFitness, fitness_child1)) {
					bestFitness = fitness_child1; // Copy best fitness, index, and position vector
					bestParticleIndex = toMutate.get(i);
					if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
					particles[bestParticleIndex].copyPosition(bestPosition);
				}
			}
			else if (fitnessFunction.isBetterThan(parent_fitness, fitness_child2) && fitnessFunction.isBetterThan(fitness_child1, fitness_child2)) {
				particles[toMutate.get(i)].position = pos_child2;

				if (fitnessFunction.isBetterThan(bestFitness, fitness_child2)) {
					bestFitness = fitness_child2; // Copy best fitness, index, and position vector
					bestParticleIndex = toMutate.get(i);
					if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
					particles[bestParticleIndex].copyPosition(bestPosition);
				}
			}
			
		}
		
	}
	//3 offsprings
	private void crossover2(ArrayList<Integer> toMutate, ArrayList<Integer> notToMutate) {

		int dimen = particles[0].position.length;
		int mut_size = notToMutate.size();
		
		for(int i=0; i<toMutate.size(); i++) {
			int rand = (int)Math.random()*mut_size;
			
			//parents' position
			double[] pos_parent1 = particles[notToMutate.get(rand)].position;		//better parent
			double[] pos_parent2 = particles[toMutate.get(i)].position;				//worse  parent
			
			//offsprings' position
			double[] pos_child1 = new double[dimen];
			double[] pos_child2 = new double[dimen];
			double[] pos_child3 = new double[dimen];
			
			int part = dimen/3;
			for(int j=0; j<dimen; j++) {
				pos_child1[j] = pos_parent1[j];
				pos_child2[j] = pos_parent1[j];
				pos_child3[j] = pos_parent1[j];
			}
			
			for(int j=0; j<part; j++) {
				pos_child1[j] = pos_parent2[j];
			}
			for(int j=part; j<2*part; j++) {
				pos_child2[j] = pos_parent2[j];
			}
			for(int j=2*part; j<dimen; j++) {
				pos_child3[j] = pos_parent2[j];
			}
			

			double fitness_child1 = fitnessFunction.evaluate(pos_child1);
			double fitness_child2 = fitnessFunction.evaluate(pos_child2);
			double fitness_child3 = fitnessFunction.evaluate(pos_child3);
			double parent_fitness = allFitness[toMutate.get(i)];

			if (fitnessFunction.isBetterThan(parent_fitness, fitness_child1) && fitnessFunction.isBetterThan(fitness_child2, fitness_child1) && fitnessFunction.isBetterThan(fitness_child3, fitness_child1)) {
				particles[toMutate.get(i)].position = pos_child1;

				if (fitnessFunction.isBetterThan(bestFitness, fitness_child1)) {
					bestFitness = fitness_child1; // Copy best fitness, index, and position vector
					bestParticleIndex = toMutate.get(i);
					if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
					particles[bestParticleIndex].copyPosition(bestPosition);
				}
			}
			else if (fitnessFunction.isBetterThan(parent_fitness, fitness_child2) && fitnessFunction.isBetterThan(fitness_child1, fitness_child2) && fitnessFunction.isBetterThan(fitness_child3, fitness_child2)) {
				particles[toMutate.get(i)].position = pos_child2;

				if (fitnessFunction.isBetterThan(bestFitness, fitness_child2)) {
					bestFitness = fitness_child2; // Copy best fitness, index, and position vector
					bestParticleIndex = toMutate.get(i);
					if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
					particles[bestParticleIndex].copyPosition(bestPosition);
				}
			}
			else if (fitnessFunction.isBetterThan(parent_fitness, fitness_child3) && fitnessFunction.isBetterThan(fitness_child1, fitness_child3) && fitnessFunction.isBetterThan(fitness_child2, fitness_child3)) {
				particles[toMutate.get(i)].position = pos_child3;

				if (fitnessFunction.isBetterThan(bestFitness, fitness_child3)) {
					bestFitness = fitness_child3; // Copy best fitness, index, and position vector
					bestParticleIndex = toMutate.get(i);
					if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
					particles[bestParticleIndex].copyPosition(bestPosition);
				}
			}
			
		}
		
	}

	void mutate(ArrayList<Integer> toMutate, ArrayList<Integer> notToMutate) {
		int dimen = particles[0].position.length;
		double[] maxPosition = new double[dimen];
		double[] minPosition = new double[dimen];

		System.out.println(toMutate.size() + ",  " + notToMutate.size());
		
		for(int j=0; j<dimen; j++) {
			maxPosition[j] = 0;
			minPosition[j] = 7; 
		}
		
		for(int i=0; i<notToMutate.size(); i++) {
			int pos = notToMutate.get(i);
			for(int j=0; j<dimen; j++) {
				double tmp = particles[pos].position[j]; 
				if(tmp+1>maxPosition[j] && tmp+1!=8)
					maxPosition[j] = (int)(tmp+1);
				else if(tmp-1<minPosition[j] && tmp-1>=0)
					minPosition[j] = (int)(tmp-1);
			}
		}
		for(int i=0; i<dimen; i++) {
			System.out.println(minPosition[i] + " - " + maxPosition[i]);
		}
		
		for(int i=0; i<toMutate.size(); i++) {
			int pos = toMutate.get(i);
			for(int j=0; j<dimen; j++) {
				double tmp = particles[pos].position[j]; 
				if(tmp>maxPosition[j] || tmp<minPosition[j]) {
					particles[pos].position[j] = (maxPosition[j] - minPosition[j]) * Math.random() + minPosition[j];
				}
			}
			
			double fit = fitnessFunction.evaluate(particles[pos]);
			allFitness[pos] = fit;
			
			if (fitnessFunction.isBetterThan(bestFitness, fit)) {
				bestFitness = fit; // Copy best fitness, index, and position vector
				bestParticleIndex = pos;
				if (bestPosition == null) bestPosition = new double[sampleParticle.getDimension()];
				particles[bestParticleIndex].copyPosition(bestPosition);
			}
		}
		
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
		double totalFitness = averageFitness*N;
		int dimension = particles[0].getPosition().length;
		
		double worstFitness=-1000000000,bestFitness=1000000000;
		for(int i=0;i < N; i++){
			System.out.println(particles[i].getFitness() + ", " + allFitness[i]);
//			bestFitness=Math.max(bestFitness, particles[i].getFitness());
//			worstFitness=Math.min(worstFitness, particles[i].getFitness());
			bestFitness=Math.min(bestFitness, allFitness[i]);
			worstFitness=Math.max(worstFitness, allFitness[i]);
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
			imass[i]=(totalFitness - allFitness[i])/(totalFitness - bestFitness);
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
//				   double temp = 0.5 * ( getG(it,mx_it) *gmass[i]*gmass[j]* direction(pos2[k], pos[k]) / euclid);
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

	   SwarmGSA_GA.accelerationArray = acc;
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
		double alfa=20,g0=1;
//		double G=g0*(Math.exp(-1*alfa*it/mx_it));
		double G=g0*(1-it/(alfa*mx_it));
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
	public double[][] particle_position() {
		double x[][] =  {{2.057666565107286,4.118057523161935,3.7946575601332886,2.606866206380526,6.968308857779919,6.013689621847535,3.279538112747324,4.916787090708778,3.2788355146494434,1.4029325902103498},
				{0.8700547532565995,1.01803892328854,6.021267077643298,5.5678072967047925,1.077891650892501,2.5209193991532186,5.629516909842428,1.0251090180204019,5.982632024016022,3.2077781176709133},
				{2.728910895733335,6.156837070933122,1.3865891429189832,6.220189761103506,1.016582110950531,1.0267067647707813,3.882948136047533,6.217413695615653,6.6185771662776745,5.629663148940408},
				{5.872952091364573,4.742803746961395,4.632748107412832,6.2617446223315,6.312989401251139,1.7721346365304358,6.4834992244963106,1.4818210591040786,1.8113000336980616,4.99887121056346},
				{4.225197499126636,6.264665910589738,4.6266395809630705,3.8362059453053083,3.0355514171431253,0.41649734632063773,4.147605267296672,1.3842486649679446,5.838500482461729,3.1395542930770883},
				{3.55144760335322,3.723055677408604,0.2881546087879775,6.791454838367747,3.5550533739908774,3.84823872549049,4.23106038817459,5.388010312515333,1.4784106870273326,5.365702485590262},
				{5.243566518825816,1.244852220188126,5.242226729962116,1.7056112404462125,1.9143019453817889,3.3436559331147544,1.4164897567336148,1.1029894690779654,3.20877169546918,6.542925737508495},
				{2.0256334282623873,0.9234189964748164,6.5138185512012425,6.631235862849848,5.498579174248739,2.087096599332509,3.814191174463874,3.4425021289812743,6.928296870040878,1.3630730882928057},
				{0.9016992326883964,1.818285498095611,3.675646888031957,1.16779797263441,2.5594320006527944,0.11131378966018646,6.316126104813711,6.208149222973486,1.2992544491495437,5.1773329599976545},
				{3.9857653036381455,6.858933086478379,5.460292748572871,0.3990806232857732,4.333786400341176,3.558391665374802,5.299356779747642,2.2266050128052486,6.534750386638272,1.742561988678084},
				{0.12124558528856,3.002330994351669,2.1517875668993645,1.9203858340453905,4.507618930188876,1.4672045943108776,0.07512862203525028,2.4162716707656013,3.307613671504124,4.31592879297684},
				{2.286509255740267,0.0753668279235995,0.2516247383268757,4.484729800779469,5.8201512716240265,0.8391035024338945,1.8416611584520832,0.5142764680054428,4.004882060378875,4.55463244191883},
				{4.557433369848086,4.333302421589327,5.879572267577871,2.681569756609567,4.39681619337608,0.15516109351359064,6.575444005579749,6.70302998569846,3.3425759825056227,6.869974457808208},
				{1.4539110091019316,2.0352283439013954,3.5142449734164805,4.349213689223148,3.8774883842081107,1.0609864546969223,0.5247615557428104,1.4574824410106757,6.6739172870551,6.587284066782197},
				{5.9159953310166005,4.287950476777796,4.630025708215806,2.5455487078129337,4.495512363717207,5.203317204249212,2.7219579866052954,5.115686294560519,3.93892909020135,5.673363528216114},
				{2.2522749516162395,0.660904608971165,2.1579486570398463,5.653078556406786,4.622043992362674,1.9101900656480812,0.6636615203839868,0.34203152260542,2.117404922018288,0.7188520477846746},
				{2.4779377636800346,6.606083294154295,2.8352150163554857,6.474851091874177,0.06272738782309728,3.126717597784531,5.71124370273863,2.799686640032211,2.2478308772889606,5.471195533354024},
				{5.547003641766627,1.6174215148907707,1.2421911490745279,0.7484236050194479,6.957675489594681,1.9339741365676053,1.5120022740058945,5.330847429054094,3.2836944774033663,5.781118088962966},
				{2.6009213602551062,4.281870413377415,1.722442712243521,3.735678915070075,1.8534127135323817,5.631220361150357,1.2440987248910718,4.224227397563827,2.481746585692534,6.487371433306716},
				{6.257072503398636,2.4671160957205087,6.490413903854821,1.1665539768579003,0.4835629220369334,1.4030214991523522,5.638713766397391,4.0520695799223665,5.094647544106522,3.2661955781227796},
				{0.8830131255435656,6.775552593373932,6.836930338962283,0.6634484530913756,4.594664506539395,1.5335053148716435,2.659137479941448,4.016427574563651,6.297818389331246,0.9178573428013808},
				{1.2171734253710442,6.097070130156989,4.713351750831704,6.548485100660509,4.6000119272276265,1.0901587192412898,6.530880145808654,5.048553197636567,2.6764165678638765,3.1667127495234726},
				{2.560312780391985,6.5002195230908795,1.6138527766149624,0.9455057730068023,4.737649973094954,6.458679415352467,0.36606771212774114,2.457301764364364,3.876718235818645,2.8866977558325098},
				{2.264685404957888,0.30745793294918544,0.06927544271849273,0.9997122796105922,1.2827780567563978,1.6573483633459039,1.2404167979894931,5.983270178328585,2.064745303279985,4.413614505597792},
				{4.56444830385889,6.859987790013934,1.4369468610636804,2.6265561330637,1.2383921582314523,6.4424052473641,1.8274024959393054,0.5550281707617237,5.496026072805817,0.3110357334831164}};
		return x;
	}public double[][] particle_velocity() {
		double x[][] =  {	{0.029658107779525822,0.023315244306486307,0.03930246053240938,0.032585588011917244,-0.011126774026204614,-0.013007303480419646,0.06970041723815701,-0.037881402937669666,0.09587644429364356,-0.052271687805914624},
				{0.07048535949395765,0.06866211857266025,-0.04477193796472992,0.07466281593584359,-0.07932186912965822,0.02960647968883673,0.07549751691417733,0.04748205665633373,-0.04989729882967753,-0.04238390143841591},
				{0.08744553386415263,0.049321678790319484,0.06782918409581581,0.03874752959813921,-0.03378658147870124,0.03872319740835778,-0.03740419881274444,-0.005529351322129658,-0.04476933064194377,-0.09905949160799937},
				{-0.027762631094139484,0.046114628143114594,0.05224742225774842,0.06190382706417891,-0.06503175521458812,-0.07860673235956017,0.03684352498701582,-0.06259439618918441,-0.09256210153664667,0.010047037702272263},
				{0.09396137978041441,-0.014418950914899503,0.06279769895485698,0.04927197753370571,-0.08054507679456975,-0.0674864647631081,-0.033412319256868536,0.09948847707291444,-0.018657423646083388,-0.05119699165064966},
				{-0.030927504374204773,0.05713566502247411,-0.045602835273849895,-0.03080661660682664,0.049725977081809314,-0.056439933762112565,0.04590117363511406,-0.032137358019644,-0.09711587038331149,-0.09230008130751358},
				{-0.013158931094540291,0.028737272035143813,-0.036843414759599666,0.08514937705995493,0.0755095982401621,0.08320032980447764,0.035006973738302055,-0.03268036330393423,0.029430481656039092,0.026820458217629767},
				{0.012982892964175693,-0.04820313405410037,-0.07328498560638197,-0.027763277098255504,0.057396475454252305,0.06862957531925964,-5.944241755511948E-4,0.04709222191185408,-0.0751730062028255,-0.09283908837666005},
				{0.07905626110296346,0.0039170056163415345,-0.0563629861830445,-0.09584582327506566,0.08874618924846744,-0.026881199718796725,-0.02823583505780118,-0.0036573091409995706,-0.07408361529253911,-0.03925217968344002},
				{-0.022629239944293375,-0.010832146061206457,-0.01313361911043838,-0.030540354548216575,0.04250387359173352,0.048352209138737334,0.04644206528832709,-0.0958918435354417,0.005738727986367723,0.08208701917871264},
				{0.02539156201452683,-0.08734181298876553,0.031848675595617826,0.04899667294446003,-0.046938760847506124,-0.09138465463120336,-0.07513474817884942,6.751453530701446E-4,-0.03545063356969831,-0.0770715518926646},
				{-0.01453785542129779,-0.052527531505816044,0.05913003855312088,0.041327553076573986,0.05441723401016571,0.022341990366275127,-0.009294456269695744,-0.002681317031175162,0.08097717004234328,-0.09005797552256044},
				{0.06124044175621429,0.011563893408021747,-0.09173072000618324,-0.09318164100620883,-0.05542387252685599,-0.004866688295720234,-0.022524176342268512,0.08337675957728002,-0.019280214091045675,-0.032876486644152034},
				{-0.08169389726915258,0.02948760498054498,0.049647122018922796,0.07639442897995388,-0.08462795861217923,-0.017052567455752926,-0.04114657593962301,-0.03731518095421481,-0.08160666376713813,-0.03634489939647505},
				{0.060332949513165396,-0.055875913563504125,-0.0012715675954358685,0.06624305836063538,-0.049638197530349815,0.07616613338727574,0.04554808285559711,-0.09996305361208487,-5.9505902327113E-4,0.0912286418367266},
				{-0.05796164336002041,-0.0924002693878444,-0.0903368632654221,-0.0011709376121763065,-0.034659590426309506,-0.05825081452952703,-0.0957789361821016,-0.09723358794270157,0.01350464303251156,0.035820535957079014},
				{-0.08554696708896986,-0.025094277320251365,-0.021442943105540263,-0.03692125685456796,0.0512881861162969,-0.08562013519811305,-0.04315388436081402,-0.051432175185958086,-0.07886778125511294,0.01653122197171142},
				{0.01615858606825138,-0.02968264838371025,-0.021311874870887662,-0.022326385539055166,-0.07633681506391672,0.06670911207544325,0.09661367605289861,-0.030237501219852844,-0.020438746417150627,-0.04706027033740659},
				{0.07589868351172213,-0.06668051767235819,0.09839380409796183,0.041897837694942874,-0.07895956078817179,0.05746487888030544,-0.06672514798945531,0.06340723403171675,-0.03945901820202063,0.04266984040142108},
				{-0.05210901053497255,0.08067023264378842,-0.09405475798619496,-0.02515246023134872,-0.09307228470172567,0.013878259715321334,-0.08841571085787553,0.07866050047891202,0.07765315581384993,-0.0016604456428520181},
				{0.022383368744336815,0.06342915861994067,-0.08628977259278939,0.028106616351307567,0.07853964447630957,0.026784434818855007,-0.07326141366017272,-0.034031076245868194,-0.01643014551199569,-0.06078502175590892},
				{0.04226800706023495,-0.06319084953963101,0.0361734864769761,-0.003082482934911715,-0.016604835711438753,-0.0362959496955773,-0.023447670576470903,8.665650360454719E-4,-0.07092393272945774,0.05478387205039348},
				{-0.09572745573906133,-0.0032224198549341693,0.05343675915923429,0.0037470608339197076,0.08745080261496349,0.09364545617084663,-0.03167482245589388,-0.07194021473763285,-0.09800872872059237,0.06779801410062636},
				{-0.036602618993693145,0.036252651487090615,-0.08520445553987276,0.054912040072577856,0.06635167527518665,0.08935973076485632,0.05631301919650772,-0.06580281880302671,-0.06950554149027816,-0.09857206499937039},
				{-0.08450562924637889,-0.08303919368537326,-0.0023580144799497937,0.010157473843394577,0.03415626551296386,-0.033975148104923816,0.08143243033585723,-0.058882938854915225,0.07796060226448229,-0.007026692510189017}};
		return x;
	}
	
}
//450053.85111111094
//36388.33
//TODO cloudlet parameters
//cloudletLength
//cloudletFileSize
//numberOfPes
