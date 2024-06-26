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

package de.unibonn.realkd.algorithms.sampling;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.function.Function;

import com.google.common.collect.ForwardingSortedSet;

import de.unibonn.realkd.algorithms.AbstractMiningAlgorithm;
import de.unibonn.realkd.algorithms.AlgorithmCategory;
import de.unibonn.realkd.algorithms.common.PatternOptimizationFunction;
import de.unibonn.realkd.data.propositions.Proposition;
import de.unibonn.realkd.data.propositions.PropositionalContext;
import de.unibonn.realkd.patterns.Pattern;
import de.unibonn.realkd.patterns.logical.LogicalDescriptor;
import de.unibonn.realkd.patterns.logical.LogicalDescriptors;
import edu.uab.consapt.sampling.TwoStepPatternSampler;
import ua.ac.be.mime.plain.PlainItem;
import ua.ac.be.mime.plain.PlainItemSet;

/**
 * Algorithm that wraps consapt based sampling of patterns with subsequent local
 * optimization (via pruning).
 * 
 * @author Mario Boley
 * @author Sandy Moens
 * @author Pavel Tokmakov
 * 
 * @version 0.6.0
 * 
 * @since 0.1.0
 *
 * @param <T>
 *            the type of patterns produced by this algorithm
 */
public class ConsaptBasedSamplingMiner<T extends Pattern<?>> extends AbstractMiningAlgorithm<T> {

	/**
	 * Mutable result data structure with limited size that allows efficient
	 * elimination of duplicates as well as efficient expulsion of worst element
	 * (in order to update with better element).
	 *
	 */
	private static class ResultSet<T extends Pattern<?>> extends ForwardingSortedSet<T> {

		private final int maximumSize;

		private final TreeSet<T> delegate;

		public ResultSet(Comparator<Pattern<?>> order, int maximumSize) {
			/*
			 * have to wrap comparator to break ties of non-identical objects
			 * because otherwise TreeSet would discard non-identical patterns
			 * with same quality
			 */
			Comparator<? super T> identityAwareComparator = (p, q) -> {
				int primaryComparison = order.compare(p, q);
				return (primaryComparison != 0) ? primaryComparison : p.equals(q) ? 0 : p.hashCode() - q.hashCode();
			};
			this.delegate = new TreeSet<>(identityAwareComparator);
			this.maximumSize = maximumSize;
		}

		@Override
		protected SortedSet<T> delegate() {
			return delegate;
		}

		@Override
		public boolean add(T p) {
			if (size() < maximumSize) {
				return delegate.add(p);
			}
			T last = delegate.last();
			if (delegate.comparator().compare(last, p) <= 0) {
				return false;
			}
			delegate.pollLast();
			return delegate.add(p);
		}

	}

	// private static final Logger LOGGER =
	// Logger.getLogger(ConsaptBasedSamplingMiner.class.getName());

	private final Function<LogicalDescriptor, ? extends T> toPattern;
	private final Function<? super T, LogicalDescriptor> toDescriptor;
	private final Integer numberOfResults;
	private final Integer numberOfSeeds;
	private final SinglePatternPostProcessor postProcessor;
	private final PropositionalContext propositionalLogic;
	private final PatternOptimizationFunction targetFunction;
	private final TwoStepPatternSampler sampler;
	private final BiFunction<LogicalDescriptor, ? super T, ? extends T> toPatternWithPrevious;

	public ConsaptBasedSamplingMiner(Function<LogicalDescriptor, T> toPattern,
			Function<Pattern<?>, LogicalDescriptor> toDescriptor, PropositionalContext propositionalLogic,
			TwoStepPatternSampler consaptSampler, PatternOptimizationFunction targetFunction,
			SinglePatternPostProcessor postProcessor, Integer numberOfResults, Integer numberOfSeeds) {
		this(toPattern, (d, p) -> toPattern.apply(d), toDescriptor, propositionalLogic, consaptSampler, targetFunction,
				postProcessor, numberOfResults, numberOfSeeds);
	}

	public ConsaptBasedSamplingMiner(Function<LogicalDescriptor, T> toPattern,
			BiFunction<LogicalDescriptor, ? super T, ? extends T> toPatternWithPrevious,
			Function<? super T, LogicalDescriptor> toDescriptor, PropositionalContext propositionalLogic,
			TwoStepPatternSampler consaptSampler, PatternOptimizationFunction targetFunction,
			SinglePatternPostProcessor postProcessor, Integer numberOfResults, Integer numberOfSeeds) {

		this.toPattern = toPattern;
		this.toPatternWithPrevious = toPatternWithPrevious;
		this.toDescriptor = toDescriptor;

		this.propositionalLogic = propositionalLogic;
		this.postProcessor = postProcessor;
		this.sampler = consaptSampler;
		this.targetFunction = targetFunction;
		this.numberOfResults = numberOfResults;
		this.numberOfSeeds = numberOfSeeds;
	}

	public final PropositionalContext getPropositionalLogic() {
		return propositionalLogic;
	}

	@Override
	protected void onStopRequest() {
		if (this.sampler != null) {
			this.sampler.setStop(true);
		}
	}

	public int numberOfResults() {
		return numberOfResults;
	}

	public final Collection<T> concreteCall() {
		ResultSet<T> results = new ResultSet<T>(targetFunction.preferenceOrder(), numberOfResults);
		int seedCounter = 0;
		int numSeeds = numberOfSeeds;

		while (!stopRequested() && seedCounter < numSeeds) {

			// sample plain set of proposition indices through unsafe API
			PlainItemSet plainPattern = sampler.getNext();
			if (plainPattern == null || stopRequested()) {
				continue;
			}

			List<Proposition> propositions = new ArrayList<>();
			for (PlainItem item : plainPattern) {
				Proposition proposition = propositionalLogic.propositions().get(Integer.parseInt(item.getName()));
				propositions.add(proposition);
			}

			LogicalDescriptor description = LogicalDescriptors.create(propositionalLogic.population(), propositions);
			T pattern = toPattern.apply(description);
			T pruned = postProcessor.prune(pattern, targetFunction, toPatternWithPrevious, toDescriptor);
			results.add(pruned);
			seedCounter++;
		}

		return results;
	}

	@Override
	public String caption() {
		return "Consapt based sampling algorithm";
	}

	@Override
	public String description() {
		return "General consapt based sampling algorithm that can be wrapped by more specific implementations";
	}

	@Override
	public AlgorithmCategory getCategory() {
		return AlgorithmCategory.OTHER;
	}

}