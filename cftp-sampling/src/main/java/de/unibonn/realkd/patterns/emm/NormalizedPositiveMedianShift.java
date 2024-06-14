/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014-2017 The Contributors of the realKD Project
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
 *
 */
package de.unibonn.realkd.patterns.emm;

import static de.unibonn.realkd.common.base.Identifier.id;

import com.fasterxml.jackson.annotation.JsonCreator;

import de.unibonn.realkd.common.KdonTypeName;
import de.unibonn.realkd.common.base.Identifier;
import de.unibonn.realkd.data.table.attribute.MetricAttribute;
import de.unibonn.realkd.patterns.models.mean.MetricEmpiricalDistribution;

/**
 * @author Mario Boley
 * 
 * @since 0.6.0
 * 
 * @version 0.6.0
 *
 */
@KdonTypeName("posNormMedianShift")
public enum NormalizedPositiveMedianShift implements UnivariateMetricEmpiricalModelDeviationMeasure {

	NORMALIZED_POSITIVE_MEDIAN_SHIFT;

	@JsonCreator
	private NormalizedPositiveMedianShift instance() {
		return NORMALIZED_POSITIVE_MEDIAN_SHIFT;
	}
	
	@Override
	public String caption() {
		return "norm. pos. median shift";
	}

	@Override
	public String description() {
		return "Difference between median target attribute value in subgroup and median target attribute value in global population (divided by global max value minus global median value) or zero of difference is negative.";
	}

	@Override
	public Identifier identifier() {
		return id("normalized_positive_median_shift");
	}

	@Override
	public double value(MetricEmpiricalDistribution refModel, MetricEmpiricalDistribution localModel,
			MetricAttribute attribute) {
		double distance = Math.max(0, localModel.medians().get(0) - refModel.medians().get(0));
		return distance / (attribute.max() - attribute.median());
	}

	@Override
	public String toString() {
		return caption();
	}

}