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

package de.unibonn.realkd.algorithms.emm;

import static com.google.common.collect.ImmutableList.copyOf;
import static de.unibonn.realkd.common.base.Identifier.id;
import static de.unibonn.realkd.patterns.emm.ExceptionalModelMining.emmPattern;
import static de.unibonn.realkd.patterns.emm.ExceptionalModelMining.extensionDescriptorToControlledEmmPatternMap;
import static de.unibonn.realkd.patterns.emm.ExceptionalModelMining.extensionDescriptorToEmmPatternMap;
import static de.unibonn.realkd.patterns.subgroups.Subgroups.controlledSubgroup;

import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;

import com.google.common.collect.ImmutableList;

import de.unibonn.realkd.algorithms.AbstractMiningAlgorithm;
import de.unibonn.realkd.algorithms.AlgorithmCategory;
import de.unibonn.realkd.algorithms.common.MiningParameters;
import de.unibonn.realkd.algorithms.common.PatternOptimizationFunction;
import de.unibonn.realkd.algorithms.sampling.ConsaptBasedSamplingMiner;
import de.unibonn.realkd.algorithms.sampling.DiscriminativityDistributionFactory;
import de.unibonn.realkd.algorithms.sampling.DistributionFactory;
import de.unibonn.realkd.algorithms.sampling.FrequencyDistributionFactory;
import de.unibonn.realkd.algorithms.sampling.SamplingParameters;
import de.unibonn.realkd.algorithms.sampling.SinglePatternPostProcessor;
import de.unibonn.realkd.algorithms.sampling.WeightedDiscriminativityDistributionFactory;
import de.unibonn.realkd.common.base.ValidationException;
import de.unibonn.realkd.common.parameter.Parameter;
import de.unibonn.realkd.common.parameter.Parameters;
import de.unibonn.realkd.common.parameter.RangeEnumerableParameter;
import de.unibonn.realkd.common.workspace.Workspace;
import de.unibonn.realkd.data.propositions.AttributeBasedProposition;
import de.unibonn.realkd.data.propositions.Proposition;
import de.unibonn.realkd.data.propositions.PropositionalContext;
import de.unibonn.realkd.data.table.DataTable;
import de.unibonn.realkd.data.table.attribute.Attribute;
import de.unibonn.realkd.data.table.attribute.MetricAttribute;
import de.unibonn.realkd.patterns.emm.ExceptionalModelMining;
import de.unibonn.realkd.patterns.emm.ExceptionalModelPattern;
import de.unibonn.realkd.patterns.emm.ModelDeviationMeasure;
import de.unibonn.realkd.patterns.logical.LogicalDescriptor;
import de.unibonn.realkd.patterns.models.Model;
import de.unibonn.realkd.patterns.models.ModelFactory;
import de.unibonn.realkd.patterns.models.mean.MetricEmpiricalDistributionFactory;
import de.unibonn.realkd.patterns.models.table.ContingencyTableModelFactory;
import de.unibonn.realkd.patterns.subgroups.ControlledSubgroup;
import de.unibonn.realkd.patterns.subgroups.ReferenceDescriptor;
import de.unibonn.realkd.patterns.subgroups.Subgroups;

/**
 * Wraps a Consapt sampling miner for randomized exceptional subgroup discovery.
 * 
 * @author Mario Boley
 * 
 * @author Sandy Moens
 * 
 * @since 0.1.0
 * 
 * @version 0.3.0
 *
 */
public final class ExceptionalSubgroupSampler extends AbstractMiningAlgorithm<ExceptionalModelPattern> {

	private final Parameter<PropositionalContext> propositionalLogicParameter;
	private final Parameter<List<Attribute<?>>> targets;
	private final Parameter<Optional<? extends Attribute<?>>> controlVariable;
	private final ModelClassParameter modelClassParameter;
	private final RangeEnumerableParameter<ModelDeviationMeasure> distanceFunctionParameter;
	private final Parameter<DataTable> dataTableParameter;
	private final Parameter<PatternOptimizationFunction> targetFunctionParameter;
	private final Parameter<Set<Attribute<?>>> descriptorAttributesParameter;
	private final Parameter<DistributionFactory> distributionFactoryParameter;
	private final Parameter<SinglePatternPostProcessor> postProcessorParameter;
	private final Parameter<Integer> numberOfResultsParameter;
	private final Parameter<Integer> numberOfSeedsParameter;

	private ConsaptBasedSamplingMiner<ExceptionalModelPattern> sampler = null;

	private final List<DistributionFactory> distributionOptions;
	private final FrequencyDistributionFactory frequencyOption;
	private final DiscriminativityDistributionFactory discriminativityOption;
	// private final DistributionOption freqTimesDiscriminativityOption;
	// private final DistributionOption posFreqSquareTimesNegInvFreqOption;
	// private final DistributionOption discriminativitySquareOption;
	private final WeightedDiscriminativityDistributionFactory weightedDiscriminativityOption;

	private ExceptionalSubgroupSampler(Workspace workspace) {
		this.dataTableParameter = MiningParameters.dataTableParameter(workspace);
		this.targets = EMMParameters.getEMMTargetAttributesParameter(dataTableParameter);
		this.modelClassParameter = new ModelClassParameter(targets);
		this.distanceFunctionParameter = EMMParameters.distanceFunctionParameter(modelClassParameter);
		this.controlVariable = EMMParameters.controlAttributeParameter(dataTableParameter, targets);
		this.descriptorAttributesParameter = EMMParameters.getEMMDescriptorAttributesParameter(dataTableParameter,
				targets, controlVariable);
		this.propositionalLogicParameter = MiningParameters.matchingPropositionalLogicParameter(workspace,
				dataTableParameter);

		// TODO:
		// filter out object if at least one target
		// value is missing
		Predicate<Proposition> targetFilter = new TargetAttributePropositionFilter(dataTableParameter, targets);
		Predicate<Proposition> controlFilter = prop -> !((prop instanceof AttributeBasedProposition) && controlVariable
				.current().map(c -> c == ((AttributeBasedProposition<?>) prop).attribute()).orElse(false));
		Predicate<Proposition> additionalPropFilter = prop -> !((prop instanceof AttributeBasedProposition)
				&& descriptorAttributesParameter.current().contains(((AttributeBasedProposition<?>) prop).attribute()));

		Predicate<Proposition> propositionFilter = additionalPropFilter.and(controlFilter).and(targetFilter);

		frequencyOption = new FrequencyDistributionFactory(propositionFilter);
		discriminativityOption = new DiscriminativityDistributionFactory(dataTableParameter, targets,
				propositionFilter);
		// freqTimesDiscriminativityOption = option(new
		// DiscriminativityDistributionFactory(dataTableParameter, targets,
		// 1, 1, 1, propositionFilter));
		// posFreqSquareTimesNegInvFreqOption = option(new
		// DiscriminativityDistributionFactory(dataTableParameter, targets,
		// 0, 2, 1, propositionFilter));
		// discriminativitySquareOption = option(new
		// DiscriminativityDistributionFactory(dataTableParameter, targets,
		// 0, 2, 2, propositionFilter));
		weightedDiscriminativityOption = new WeightedDiscriminativityDistributionFactory(dataTableParameter, targets,
				propositionalLogicParameter, 0, 2, 1, new RowWeightComputer.UniformRowWeightComputer(),
				propositionFilter);
		this.distributionOptions = ImmutableList.of(frequencyOption, discriminativityOption,
				weightedDiscriminativityOption);

		this.distributionFactoryParameter = Parameters.rangeEnumerableParameter(id("seed_dist"), "Seed distribution",
				"The probability distribution on the pattern space that is used to generate random seeds for EMM pattern search",
				DistributionFactory.class, () -> this.distributionOptions, targets, propositionalLogicParameter,
				descriptorAttributesParameter);
		this.targetFunctionParameter = EMMParameters.emmTargetFunctionParameter(modelClassParameter, controlVariable);
		this.postProcessorParameter = SamplingParameters.postProcessingParameter();
		this.numberOfResultsParameter = SamplingParameters.numberOfResultsParameter();
		this.numberOfSeedsParameter = SamplingParameters.numberOfSeedsParameter(numberOfResultsParameter);

		registerParameter(dataTableParameter);
		registerParameter(targets);
		registerParameter(modelClassParameter);
		registerParameter(distanceFunctionParameter);
		registerParameter(controlVariable);
		registerParameter(propositionalLogicParameter);
		registerParameter(descriptorAttributesParameter);
		registerParameter(targetFunctionParameter);
		registerParameter(numberOfResultsParameter);
		registerParameter(numberOfSeedsParameter);
		registerParameter(distributionFactoryParameter);
		registerParameter(postProcessorParameter);
	}

	// private Function<LogicalDescriptor, ExceptionalModelPattern>
	// toPatternMap() {
	// if (controlVariable.current().isPresent()) {
	// List<Attribute<?>> controlAttributes =
	// ImmutableList.of(controlVariable.current().get());
	// MetricEmpiricalDistribution referenceControlModel =
	// MetricEmpiricalDistributionFactory.INSTANCE
	// .getModel(table, controlAttributes);
	// return extensionDescriptorToControlledEmmPatternMap(table, targetAttr,
	// modelFactory, globalModel,
	// distanceMeasurementProc, ImmutableList.of(), controlAttributes,
	// referenceControlModel,
	// MetricEmpiricalDistributionFactory.INSTANCE);
	// } else {
	// return extensionDescriptorToEmmPatternMap(table, targetAttr,
	// modelFactory, globalModel,
	// distanceMeasurementProc, ImmutableList.of());
	// }
	// }

	private static class EmmPatternConstructionSetting<M extends Model> {

		public final DataTable table;
		public final List<Attribute<?>> targetAttributes;
		public final ModelFactory<? extends M> modelFactory;
		public final ModelDeviationMeasure deviationMeasure;

		private Model referenceTargetModel = null;

		public EmmPatternConstructionSetting(DataTable table, List<Attribute<?>> targetAttributes,
				ModelFactory<? extends M> modelFactory, ModelDeviationMeasure distanceMeasurementProc) {
			this.table = table;
			this.targetAttributes = targetAttributes;
			this.modelFactory = modelFactory;
			this.deviationMeasure = distanceMeasurementProc;
		}

		public Model referenceTargetModel() {
			if (referenceTargetModel == null) {
				referenceTargetModel = modelFactory.getModel(table, targetAttributes);
			}
			return referenceTargetModel;
		}

		public Function<LogicalDescriptor, ExceptionalModelPattern> simpleDescriptorToPatternMap() {
			return extensionDescriptorToEmmPatternMap(table, targetAttributes, modelFactory,
					d -> referenceTargetModel(), d -> ReferenceDescriptor.global(table.population()),
					deviationMeasure, ImmutableList.of());
		}

		public BiFunction<LogicalDescriptor, ExceptionalModelPattern, ExceptionalModelPattern> descriptorAndPredecessorToPatternMap() {
			return (d, p) -> {
				if (p.descriptor().extensionDescriptor().supportSet().equals(d.supportSet())) {
					return ExceptionalModelMining.emmPattern(
							Subgroups.subgroup(d, ReferenceDescriptor.global(table.population()), table,
									targetAttributes, modelFactory, referenceTargetModel, p.descriptor().localModel()),
							deviationMeasure, ImmutableList.of());
				}
				Model localModel = modelFactory.getModel(table, targetAttributes, d.supportSet());
				return ExceptionalModelMining.emmPattern(
						Subgroups.subgroup(d, ReferenceDescriptor.global(table.population()), table,
								targetAttributes, modelFactory, referenceTargetModel, localModel),
						deviationMeasure, ImmutableList.of());
			};
		}
	}

	private static class ControlVariableSetting<T extends Model, C extends Model>
			extends EmmPatternConstructionSetting<T> {

		public final List<Attribute<?>> controlAttributes;

		public final ModelFactory<C> factory;

		private C referenceControlModel = null;

		public C referenceControlModel() {
			if (referenceControlModel == null) {
				referenceControlModel = factory.getModel(table, controlAttributes);
			}
			return referenceControlModel;
		}

		@Override
		public Function<LogicalDescriptor, ExceptionalModelPattern> simpleDescriptorToPatternMap() {
			return extensionDescriptorToControlledEmmPatternMap(table, targetAttributes, modelFactory,
					referenceTargetModel(), deviationMeasure, ImmutableList.of(), controlAttributes,
					referenceControlModel(), factory);
		}

		@Override
		public BiFunction<LogicalDescriptor, ExceptionalModelPattern, ExceptionalModelPattern> descriptorAndPredecessorToPatternMap() {
			return (d, p) -> {
				if (p.descriptor().extensionDescriptor().supportSet().equals(d.supportSet())) {
					ControlledSubgroup<?, ?> controlledSubgroup = (ControlledSubgroup<?, ?>) p.descriptor();
					return emmPattern(
							controlledSubgroup(d, table, targetAttributes, modelFactory, referenceTargetModel(),
									controlledSubgroup.localModel(), controlAttributes, factory,
									controlledSubgroup.referenceControlModel(), controlledSubgroup.localControlModel()),
							deviationMeasure, ImmutableList.of());
				}
				Model localModel = modelFactory.getModel(table, targetAttributes, d.supportSet());
				Model localControlModel = factory.getModel(table, controlAttributes, d.supportSet());
				return emmPattern(
						controlledSubgroup(d, table, targetAttributes, modelFactory, referenceTargetModel(), localModel,
								controlAttributes, factory, referenceControlModel(), localControlModel),
						deviationMeasure, ImmutableList.of());
			};
		}

		public ControlVariableSetting(DataTable table, List<Attribute<?>> targetAttributes,
				ModelFactory<T> targetModelFactory, ModelDeviationMeasure targetDeviationMeasure,
				List<Attribute<?>> controlAttributes, ModelFactory<C> controlModelFactory) {
			super(table, targetAttributes, targetModelFactory, targetDeviationMeasure);
			this.controlAttributes = controlAttributes;
			this.factory = controlModelFactory;
		}

	}

	private static <M extends Model> EmmPatternConstructionSetting<M> constructionSetting(DataTable table,
			List<Attribute<?>> targetAttributes, ModelFactory<M> modelFactory,
			ModelDeviationMeasure distanceMeasurementProc, Optional<? extends Attribute<?>> controlAttribute) {
		if (controlAttribute.isPresent()) {
			List<Attribute<?>> controlAttributes = ImmutableList.of(controlAttribute.get());
			ModelFactory<?> factory = (controlAttribute.get() instanceof MetricAttribute)
					? MetricEmpiricalDistributionFactory.INSTANCE
					: ContingencyTableModelFactory.INSTANCE;
			return new ControlVariableSetting<>(table, targetAttributes, modelFactory, distanceMeasurementProc,
					controlAttributes, factory);
		} else {
			return new EmmPatternConstructionSetting<M>(table, targetAttributes, modelFactory, distanceMeasurementProc);
		}
	}

	public static ExceptionalSubgroupSampler exceptionalSubgroupSampler(Workspace workspace) {
		return new ExceptionalSubgroupSampler(workspace);
	}

	@Override
	protected Collection<ExceptionalModelPattern> concreteCall() throws ValidationException {

		final DataTable table = dataTableParameter.current();
		final List<Attribute<?>> targetAttr = targets.current();
		final ModelFactory<?> modelFactory = modelClassParameter.current().get();
		final ModelDeviationMeasure distanceMeasurementProc = distanceFunctionParameter.current();
		final Optional<? extends Attribute<?>> controlAttr = controlVariable.current();

		final EmmPatternConstructionSetting<?> setting = constructionSetting(table, targetAttr, modelFactory,
				distanceMeasurementProc, controlAttr);
		final Function<LogicalDescriptor, ExceptionalModelPattern> toPattern = setting.simpleDescriptorToPatternMap();
		final BiFunction<LogicalDescriptor, ExceptionalModelPattern, ExceptionalModelPattern> toPatternWithPrevious = setting
				.descriptorAndPredecessorToPatternMap();

		sampler = new ConsaptBasedSamplingMiner<ExceptionalModelPattern>(
				// new
				// ParameterBoundLogicalDescriptorToEmmPatternMap(dataTableParameter,
				// targets, modelClassParameter,
				// distanceFunctionParameter),
				toPattern, toPatternWithPrevious,
				pattern -> ((ExceptionalModelPattern) pattern).descriptor().extensionDescriptor(),
				propositionalLogicParameter.current(),
				distributionFactoryParameter.current().getDistribution(propositionalLogicParameter.current()),
				targetFunctionParameter.current(), postProcessorParameter.current(), numberOfResultsParameter.current(),
				numberOfSeedsParameter.current());

		Collection<ExceptionalModelPattern> result = sampler.call();
		sampler = null;
		return result;
	}

	@Override
	protected void onStopRequest() {
		if (sampler != null) {
			sampler.requestStop();
		}
	}

	public void useDiscriminativitySeedOption() {
		this.distributionFactoryParameter.set(discriminativityOption);
	}

	public DiscriminativityDistributionFactory discriminativityOption() {
		return discriminativityOption;
	}

	public void setWeightedDiscriminativityOption() {
		this.distributionFactoryParameter.set(weightedDiscriminativityOption);
	}

	public DistributionFactory getDistributionFactory() {
		return distributionFactoryParameter.current();
	}

	public List<Attribute<?>> targetAttributes() {
		return this.targets.current();
	}

	@Override
	public String toString() {
		return "EMMSampler|" + getDistributionFactory() + "|" + targetFunctionParameter.current();
	}

	@Override
	public String caption() {
		return "Randomized Exceptional Subgroup Discovery";
	}

	@Override
	public String description() {
		return "2-Step direct pattern sampling for discovering outstanding subgroups.";
	}

	@Override
	public AlgorithmCategory getCategory() {
		return AlgorithmCategory.EXCEPTIONAL_SUBGROUP_DISCOVERY;
	}

	public Parameter<PropositionalContext> propositionalLogicParameter() {
		return this.propositionalLogicParameter;
	}

	public Parameter<DataTable> dataTableParameter() {
		return dataTableParameter;
	}

	public Parameter<List<Attribute<?>>> targetAttributesParameter() {
		return targets;
	}

	public ExceptionalSubgroupSampler targetAttributes(Attribute<?>... attributes) {
		targets.set(copyOf(attributes));
		return this;
	}

	public ExceptionalSubgroupSampler targetAttributes(List<? extends Attribute<?>> attributes) {
		targets.set(copyOf(attributes));
		return this;
	}

	public ModelClassParameter modelClassParameter() {
		return modelClassParameter;
	}

	public RangeEnumerableParameter<ModelDeviationMeasure> modelDistanceFunctionParameter() {
		return distanceFunctionParameter;
	}

	public Parameter<DistributionFactory> seedDistributionParameter() {
		return this.distributionFactoryParameter;
	}

	public Parameter<Set<Attribute<?>>> descriptionAttributesParameter() {
		return descriptorAttributesParameter;
	}

	public Parameter<?> objectiveFunctionParameter() {
		return targetFunctionParameter;
	}

	public ExceptionalSubgroupSampler numberOfResults(int numberOfResults) {
		numberOfResultsParameter.set(numberOfResults);
		return this;
	}

	public Parameter<Integer> numberOfResultsParameter() {
		return numberOfResultsParameter;
	}

	public void numberOfSeeds(int i) {
		numberOfSeedsParameter.set(i);
	}

	public Parameter<Integer> numberOfSeedsParameter() {
		return numberOfSeedsParameter;
	}

	public ExceptionalSubgroupSampler useSingleEventModel() {
		modelClassParameter().set(modelClassParameter().bernoulliOption());
		return this;
	}

	public ExceptionalSubgroupSampler positiveCategory(String category) {
		modelClassParameter().bernoulliOption().positiveCategory().setByString(category);
		return this;
	}

}
