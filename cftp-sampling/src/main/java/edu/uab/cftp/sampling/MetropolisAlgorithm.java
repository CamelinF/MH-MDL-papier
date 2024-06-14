package edu.uab.cftp.sampling;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import edu.uab.cftp.sampling.distribution.OrderedBaseObject;
import edu.uab.cftp.sampling.distribution.State;
import edu.uab.cftp.sampling.distribution.StateSpaceProbabilityDistribution;

/**
 * Implements a simple Metropolis Chain simulation on a finite state space with
 * acceptance probablity p(x)/p(y), where x is the current state and y is a
 * newly proposed state and p(.) is the probability of being in a given state.
 * 
 * @author Sandy Moens
 */

public class MetropolisAlgorithm implements
		InterruptiblePerfectSamplingAlgorithm {

	protected StateSpaceProbabilityDistribution distribution;
	protected State curState = null;
	protected boolean accepted = false;

	protected double newValue = 0;
	protected double currentValue = 0;

	public MetropolisAlgorithm() {
	}

	public MetropolisAlgorithm(StateSpaceProbabilityDistribution distribution) {
		this.distribution = distribution;
		reset();
	}

	@Override
	public void reset() {
		this.curState = this.distribution.getHardJupp();
		this.currentValue = getValue(this.curState);
		this.accepted = false;
	}

	protected double getValue(State state) {
		return state.getProbability();
	}

	@Override
	public State run(long time, long randomSeed) {

		State newState;
		Random random = new Random(randomSeed);

		//System.out.println("avec la seed :"+randomSeed );

		double test;
		double valeurTest;
		boolean acc;
		//int compteur=0;
		//System.out.println((time));
		//double maxi=-1;
		//OrderedBaseObject maxState=null;
		for (int i = 0; i < time; i++) {


			newState = this.distribution.drawNextState(random);

			this.newValue = getValue(newState);
			/*if (maxi < this.newValue){
				maxi=this.newValue;
				maxState=newState.getBaseObject();
				compteur=1;
			}else{
				if(maxi == this.newValue) {
					compteur += 1;
				}
			}*/
			test =random.nextDouble();

			valeurTest= newValue / this.currentValue;

			//System.out.println(valeurTest);

			acc=false;
			if (test < valeurTest) {
				this.curState = newState;
				this.currentValue = this.newValue;
				this.accepted = true;
				acc=true;
			}
			/*if (i ==time-1 &&  time>Math.pow(2,20)){
				this.curState = newState;
				this.currentValue = this.newValue;
				this.accepted = true;
				acc=true;
			}*/
			if (acc) {
				//System.out.println(test);
				//System.out.println(" seuil de "+test+",NEWVALUE = "+ this.newValue+", CURRENT= "+this.currentValue+", Valeur du test ="+valeurTest);
				//System.out.println(" ===========> accepté");
			} else {
				//System.out.print(" seuil de "+test+",NEWVALUE = "+ this.newValue+", CURRENT= "+this.currentValue+", Valeur du test ="+valeurTest);

				//System.out.println("");
			}
		}
		//System.out.println(" valuer maxi trouvé= "+maxi+" trouvé "+compteur+" fois donc min(q/w)="+1/maxi);
		//System.out.println(maxState.toString());
		if (this.accepted) {
			//System.out.println("epoch : "+time);
			return this.curState;
		} else {
			return null;
		}
	}

	@Override
	public State run(int beg, List<long[]> randomSeeds) {
		Iterator<long[]> iterator = randomSeeds.iterator();
		for (int i = beg; i > 0; i--) {
			iterator.next();
		}
		long[] randomSeed;
		while (iterator.hasNext()) {
			randomSeed = iterator.next();
			run(randomSeed[0], randomSeed[1]);
		}

		if (this.accepted) {
			return this.curState;
		} else {
			return null;
		}
	}
}
