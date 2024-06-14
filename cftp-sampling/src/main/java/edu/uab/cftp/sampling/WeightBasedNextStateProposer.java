package edu.uab.cftp.sampling;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import edu.uab.cftp.sampling.distribution.OrderedBaseObject;
import edu.uab.cftp.sampling.distribution.State;
import ua.ac.be.mime.plain.PlainTransaction;
import ua.ac.be.mime.tool.Utils;

public class WeightBasedNextStateProposer extends NextStateProposer {

	private double[] cumulatedValues;
	private double[] values;

	@Override
	public void setTransactions(List<PlainTransaction> transactions) {
		this.transactions = transactions;
		this.cumulatedValues = new double[this.transactions.size()];
		this.values = new double[this.transactions.size()];
		computeCumulatedValues();
		/*System.out.print("Values = [");
		for(double i :this.values ){
			System.out.print(i+",");
		}
		System.out.println("]");*/
	}

	private double getWeight(PlainTransaction transaction) {
		double n =this.distribution.getWeightInDistribution(transaction.getItemsAsBitSet());

		double m=1.0 / this.distribution.cardinality();
		//System.out.println(" n = "+n+" m = "+m);
		double weight = Math.pow(n, m);
		if (Double.isNaN(weight)) {
			return 0;
		} else {

			return weight;
		}
	}

	private void computeCumulatedValues() {
		int index = 0;
		double normalizationFactor = 0;
		for (PlainTransaction transaction : this.transactions) {
			normalizationFactor += (this.values[index++] = getWeight(transaction));

		}

		double value = 0.0;
		index = 0;
		for (double v : this.values) {
			this.cumulatedValues[index++] = (value += (v / normalizationFactor));
		}
	}

	private int getIndex(double d) {
		return Utils.logIndexSearch(this.cumulatedValues, d);
	}

	@Override
	public State nextState(Random random) {
		int[] indices = new int[this.distribution.cardinality()];
		for (int index = 0; index < this.distribution.cardinality(); index++) {
			indices[index] = -1;
		}
		int x;

		double y;
		boolean encore;
		for (int index = 0; index < this.distribution.cardinality(); index++) {
			x=-1;
			encore = true;
			while (encore) {
				encore=false;
				y=random.nextDouble();
				x = getIndex(y);
				//x = random.nextInt(this.transactions.size());
				//System.out.println("avec "+y+" => "+ x);
				for (int s : indices){
					if (x == s){
						encore=true;
					}
				}

			}
			indices[index] = x;
		}
		/*System.out.print("indices=[");
		for (int index = 0; index < this.distribution.cardinality()-1; index++) {
			System.out.print(indices[index]+",");
		}
		System.out.println(indices[this.distribution.cardinality()-1]+"]");*/
		OrderedBaseObject baseObject = createBaseObject(indices);
		double w = this.distribution.baseObjectWeight(baseObject);
		//System.out.println( "\t etat proposé "+baseObject.toString()+" avec un poids de "+w);
		return new State(w,
				baseObject);
	}

	@Override
	public double getPotential(State stateFrom, State stateTo) {
		double product = 1.0;
		for (int index : stateTo.getBaseObject().getIndices()) {
			product *= this.values[index];
		}
		return product;
	}

	@Override
	public OrderedBaseObject createBaseObject(int[] indices) {
		List<PlainTransaction> transactions = new ArrayList<PlainTransaction>(
				this.distribution.cardinality());
		for (int index : indices) {
			transactions.add(this.transactions.get(index));
		}
		return new OrderedBaseObject(indices);
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
		return false;
	}
}
