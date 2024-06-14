package edu.uab.cftp.sampling;

import java.util.List;
import java.util.Random;

import edu.uab.cftp.sampling.distribution.LabeledOrderedBaseObject;
import edu.uab.cftp.sampling.distribution.OrderedBaseObject;
import edu.uab.cftp.sampling.distribution.State;
import edu.uab.cftp.sampling.distribution.old.DiscriminativityDistribution;
import edu.uab.cftp.sampling.distribution.old.FrequencyDiscriminativityDistribution;
import edu.uab.cftp.sampling.distribution.old.PosFrequencyDiscriminativityDistribution;
import ua.ac.be.mime.plain.PlainTransaction;

/**
 * 
 * @author Sandy Moens
 */

public class UniformNextStateProposer extends NextStateProposer {

	private int transactionsSize;
	private double largestValue = Double.MAX_VALUE;

	@Override
	public void setTransactions(List<PlainTransaction> transactions) {
		this.transactions = transactions;
		this.transactionsSize = this.transactions.size();
		this.largestValue = getWeight(this.transactions.get(0));
	}

	private double getWeight(PlainTransaction transaction) {
		return this.distribution.getWeightInProposal(transaction
				.getItemsAsBitSet());
	}

	@Override
	public State nextState(Random random) {
		int[] indices = new int[this.distribution.cardinality()];
		//for (int i = 0; i < this.distribution.cardinality(); i++) {
		//	indices[i] = random.nextInt(this.transactionsSize);
		//}
		for (int index = 0; index < this.distribution.cardinality(); index++) {
			indices[index] = -1;
		}
		int x;
		boolean encore;
		for (int index = 0; index < this.distribution.cardinality(); index++) {

			x=-1;
			encore = true;
			while (encore) {
				encore=false;
				x = random.nextInt(this.transactionsSize);
				for (int s : indices) {
					if (x == s){
						encore=true;
					}
				}

			}
			indices[index] = x;
		}
		OrderedBaseObject baseObject = new OrderedBaseObject(indices);


		double w = this.distribution.baseObjectWeight(baseObject);
		//System.out.print( "\t etat proposé "+baseObject.toString()+" avec un poids de "+w);
		return new State(w,
				baseObject);
	}

	@Override
	public double getPotential(State stateFrom, State stateTo) {
			return this.largestValue;
	}

	@Override
	public OrderedBaseObject createBaseObject(int[] indices) {
		if (this.distribution instanceof DiscriminativityDistribution
				|| this.distribution instanceof FrequencyDiscriminativityDistribution
				|| this.distribution instanceof PosFrequencyDiscriminativityDistribution) {
			return new LabeledOrderedBaseObject(indices);
		} else {
			return new OrderedBaseObject(indices);
		}
	}

	@Override
	public int randomBitsForNextState() {
		return this.distribution.cardinality();
	}

	@Override
	public boolean isSparse() {
		return false;
	}

	@Override
	public boolean isUniform() {
		return true;
	}
}
