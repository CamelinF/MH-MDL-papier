/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014-2018 The Contributors of the realKD Project
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
package de.unibonn.realkd.patterns.models.conditional;

import java.util.Map;
import java.util.Set;

import de.unibonn.realkd.patterns.models.Model;
import de.unibonn.realkd.patterns.models.table.Cell;

/**
 * @author Mario Boley
 * 
 * @author Kailash Budhathoki
 * 
 * @since 0.7.0
 * 
 * @version 0.7.0
 *
 */
public class DiscretelyConditionedBernoulli implements Model {

	private final Map<Cell, EmpiricalBernoulliDistribution> conditionalBernoulliTables;

	DiscretelyConditionedBernoulli(Map<Cell, EmpiricalBernoulliDistribution> conditionalBernoulliTables) {
		this.conditionalBernoulliTables = conditionalBernoulliTables;
	}

	public EmpiricalBernoulliDistribution conditionalTable(Cell cell) {
		return conditionalBernoulliTables.getOrDefault(cell, new EmpiricalBernoulliDistribution(0, 0));
	}

	public Set<Cell> cells() {
		return conditionalBernoulliTables.keySet();
	}
}
