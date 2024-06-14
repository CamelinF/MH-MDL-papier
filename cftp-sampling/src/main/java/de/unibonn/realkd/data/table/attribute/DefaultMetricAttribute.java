/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 University of Bonn
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package de.unibonn.realkd.data.table.attribute;

import static de.unibonn.realkd.common.math.types.ClosedInterval.closedInterval;
import static java.lang.Math.abs;
import static java.util.Comparator.naturalOrder;

import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;

import de.unibonn.realkd.common.IndexSet;
import de.unibonn.realkd.common.base.Identifier;
import de.unibonn.realkd.common.math.types.ClosedInterval;

public final class DefaultMetricAttribute extends DefaultAttribute<Double>
		implements OrdinalAttribute<Double>, MetricAttribute {

	private final ClosedInterval range;

	private final double mean, variance, thirdCentralMoment, avgAbsMedDev;

	private final OrderedValueContainer<Double> orderedValueContainer;

	@JsonCreator
	DefaultMetricAttribute(@JsonProperty("identifier") Identifier identifier, @JsonProperty("name") String name, @JsonProperty("description") String description,
			@JsonProperty("values") List<Double> values) {
		super(identifier, name, description, values, Double.class);
		this.orderedValueContainer = new OrderedValueContainer<>(values, naturalOrder());
		double _mean = 0.0;

		for (int i = 0; i < values.size(); i++) {
			Double d = values.get(i);
			if (d != null) {
				_mean += d;
			}
		}

		mean = _mean / (values.size() - missingPositions().size());

		double _variance = 0.0;
		double _thirdCentralMoment = 0.0;
		double _avgAbsMeanDev = 0.0;

		for (double value : nonMissingValues()) {
			_variance += (mean - value) * (mean - value);
			_thirdCentralMoment += (value - mean) * (value - mean) * (value - mean);
			_avgAbsMeanDev += Math.abs(median() - value);
		}
		int m = numberOfNonMissingValues();
		variance = (m > 0) ? _variance / m : Double.NaN;
		thirdCentralMoment = (m > 0) ? _thirdCentralMoment / m : Double.NaN;
		avgAbsMedDev = (m > 0) ? _avgAbsMeanDev / m : Double.NaN;
		range = closedInterval(min(), max());
	}

	@Override
	public Double quantile(double frac) {
		return orderedValueContainer.quantile(frac);
	}

	@Override
	public double mean() {
		return mean;
	}

	@Override
	public Double median() {
		return orderedValueContainer.median();
	}

	@Override
	public Double max() {
		return orderedValueContainer.max();
	}

	@Override
	public Double min() {
		return orderedValueContainer.min();
	}

	@Override
	public double variance() {
		return variance;
	}

	@Override
	public double averageAbsoluteMedianDeviation() {
		return avgAbsMedDev;
	}

	@Override
	public double standardDeviation() {
		return Math.sqrt(variance);
	}

	@Override
	public double skew() {
		return thirdCentralMoment / Math.pow(standardDeviation(), 3);
	}

	public List<Integer> sortedNonMissingRowIndices() {
		return orderedValueContainer.sortedNonMissingRowIndices();
	}

	public Double lowerQuartile() {
		return orderedValueContainer.lowerQuartile();
	}

	public Double upperQuartile() {
		return orderedValueContainer.upperQuartile();
	}

	@Override
	public double meanOnRows(IndexSet rowSet) {
		double result = 0.0;
		int numberOfNonMissingValues = 0;

		for (int rowIndex : rowSet) {
			if (valueMissing(rowIndex)) {
				continue;
			}
			result += value(rowIndex);
			numberOfNonMissingValues++;
		}

		if (numberOfNonMissingValues == 0) {
			return Double.NaN;
		}

		return result / numberOfNonMissingValues;
	}

	@Override
	public double averageAbsoluteMedianDeviationOnRows(IndexSet rowSet) {
		double med = meanOnRows(rowSet);
		double result = 0.0;
		int numberOfNonMissingValues = 0;

		for (int rowIndex : rowSet) {
			if (valueMissing(rowIndex)) {
				continue;
			}
			result += abs(med - value(rowIndex));
			numberOfNonMissingValues++;
		}

		if (numberOfNonMissingValues == 0) {
			return Double.NaN;
		}

		return result / numberOfNonMissingValues;
	}

	public Double medianOnRows(IndexSet rowSet) {
		return orderedValueContainer.getMedianOnRows(rowSet).orElse(Double.NaN);
	}

	@Override
	public Comparator<Double> valueComparator() {
		return orderedValueContainer.comparator();
	}

	@Override
	public int orderNumber(Double value) {
		return orderedValueContainer.orderNumber(value);
	}

	@Override
	public int orderNumberOnRows(Double value, Set<Integer> rows) {
		return orderedValueContainer.orderNumberOnRows(value, rows);
	}

	@Override
	public ClosedInterval range() {
		return range;
	}

	@Override
	public int inverseOrderNumber(Double value) {
		return orderedValueContainer.inverseOrderNumber(value);
	}
}
