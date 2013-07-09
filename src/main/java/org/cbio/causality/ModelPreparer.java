package org.cbio.causality;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.Searcher;
import org.biopax.paxtools.pattern.constraint.*;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class ModelPreparer
{
	private static final Pattern PARTICIPANT_AND_CONTROLLER_PATTERN = 
		prepareParticipantAndControllerPattern();
	
	private static final Pattern INPUT_TO_CONV_PATTERN = prepareInputToConv();
	private static final Pattern OUTPUT_TO_CONV_PATTERN = prepareOutputToConv();
	private static final Pattern BOTH_INPUT_AND_OUTPUT_PATTERN = prepareBothInputAndOutputPattern();
	private static final Pattern RELATED_ER_PATTERN = prepareRelatedERPattern();
	private static final Pattern SIMPLE_MEMBER_PATTERN = prepareSimpleMemberPattern();
	
	private static final PathAccessor COMPLEX_ACC = new PathAccessor("PhysicalEntity/componentOf*");
	private static final PathAccessor PARTICIPATES_ACC = new PathAccessor("PhysicalEntity/participantOf");
	private static final PathAccessor PARTICIPANTS_ACC = new PathAccessor("Interaction/participant");
	
	private static final Random rand = new Random();
	
	public static void prepare(Model model)
	{
		System.out.println("----Fix control and participant");
		removeControlIfAlsoParticipant(model);
		System.out.println("----Fixed");
		System.out.println("----Fix both input and output");
		fixBothInputAndOutputs(model);
		System.out.println("----Fixed");
		System.out.println("----Process loops");
		processLoops(model);
		System.out.println("----Processed");
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("/home/ozgun/Desktop/MDM2-TP53-Neigh.owl"));
		prepare(model);
		ActivityRecognizer ar = new ActivityRecognizer(model);
		ar.run();
		h.convertToOWL(model, new FileOutputStream("/home/ozgun/Desktop/MDM2-TP53-Neigh-prepared.owl"));
	}

	public static void writePatternMatches() throws FileNotFoundException
	{
		Searcher.searchInFile(BOTH_INPUT_AND_OUTPUT_PATTERN, "/home/ozgun/Desktop/cpath2.owl",
			"/home/ozgun/Desktop/pattern-matches/BOTH_INPUT_AND_OUTPUT_PATTERN.owl");
	}

	private static String generateID()
	{
		return "ModelPreparer-" + System.currentTimeMillis() + "-" + rand.nextDouble();
	}
	
	public static void removeControlIfAlsoParticipant(Model model)
	{
		for (Match match : Searcher.searchPlain(model, PARTICIPANT_AND_CONTROLLER_PATTERN))
		{
			Control c = (Control) match.get(2);
			
			if (c.getControlled().size() != 1) System.err.println("Control controls more than one. ID = " + c.getRDFId());
			if (c.getController().size() != 1) System.err.println("Control has more than one controller. ID = " + c.getRDFId());

			model.remove(c);
			for (Process controlled : new HashSet<Process>(c.getControlled()))
			{
				c.removeControlled(controlled);
			}
			for (Controller controller : new HashSet<Controller>(c.getController()))
			{
				c.removeController(controller);
			}
			for (Pathway pathway : c.getPathwayComponentOf())
			{
				pathway.removePathwayComponent(c);
			}
			for (PathwayStep step : c.getStepProcessOf())
			{
				step.removeStepProcess(c);
			}
		}
	}

	public static void fixBothInputAndOutputs(Model model)
	{
		for (Conversion cnv : new HashSet<Conversion>(model.getObjects(Conversion.class)))
		{
			for (Match match : Searcher.search(cnv, BOTH_INPUT_AND_OUTPUT_PATTERN))
			{
				PhysicalEntity pe = (PhysicalEntity) match.get(1);
				cnv.removeLeft(pe);
				cnv.removeRight(pe);
				Control ctrl = model.addNew(Control.class, generateID());
				ctrl.addController(pe);
				ctrl.addControlled(cnv);
				ctrl.setControlType(ControlType.ACTIVATION);
			}
		}
	}
	
	private static void processLoops(Model model)
	{
		Set<Set<Conversion>> groups = new HashSet<Set<Conversion>>();
		for (int i = 2; i < 5; i++)
		{
			System.out.println("running for loop size " + i);
			List<Match> matches = Searcher.searchPlain(model, prepareLoop(i));
			groups.addAll(extractUniqueConversionGroups(matches));
			System.out.println("matches.size() = " + matches.size());
			System.out.println("groups.size() = " + groups.size());
		}
		for (Set<Conversion> group : groups)
		{
			createReplacementConversion(model, group);
		}
	}
	
	private static void createReplacementConversion(Model model, Set<Conversion> cons)
	{
		Set<PhysicalEntity> inputs = new HashSet<PhysicalEntity>();
		Set<PhysicalEntity> outputs = new HashSet<PhysicalEntity>();
		Set<PhysicalEntity> effectors = new HashSet<PhysicalEntity>();

		for (Object o : PARTICIPANTS_ACC.getValueFromBeans(cons))
		{
			if (o instanceof PhysicalEntity)
			{
				PhysicalEntity pe = (PhysicalEntity) o;

				if (pe instanceof SmallMolecule) continue;

				boolean isInput = isAssociated(pe, true, cons);
				boolean isOutput = isAssociated(pe, false, cons);
				
				if (isInput && isOutput) effectors.add(pe);
				else if (isInput) inputs.add(pe);
				else if (isOutput) outputs.add(pe);
			}
		}

		if (!inputs.isEmpty() && !outputs.isEmpty() && !effectors.isEmpty())
		{
			Set<PhysicalEntity> inOut = new HashSet<PhysicalEntity>(inputs);
			inOut.addAll(outputs);

			Set<EntityReference> inOutERs = Searcher.searchAndCollect(
				inOut, RELATED_ER_PATTERN, 2, EntityReference.class);

			Conversion conv = model.addNew(Conversion.class, generateID());

			for (PhysicalEntity pe : inputs) conv.addLeft(pe);
			for (PhysicalEntity pe : outputs) conv.addRight(pe);

			Control ctrl = model.addNew(Control.class, generateID());

			for (PhysicalEntity pe : Searcher.searchAndCollect(effectors, SIMPLE_MEMBER_PATTERN, 1, PhysicalEntity.class))
			{
				if (relatedERsDifferent(pe, inOutERs)) ctrl.addController(pe);
			}

			ctrl.addControlled(conv);
		}

	}

	private static Set<Set<Conversion>> extractUniqueConversionGroups(List<Match> matches)
	{
		Set<Set<Conversion>> set = new HashSet<Set<Conversion>>();
		Set<Integer> encountered = new HashSet<Integer>();

		for (Match match : matches)
		{
			Set<Conversion> cvns = getConversions(match);
			Integer hash = hashSum(cvns);
			if (!encountered.contains(hash))
			{
				encountered.add(hash);
				set.add(cvns);
			}
		}
		return set;
	}
	
	private static Set<Conversion> getConversions(Match match)
	{
		Set<Conversion> set = new HashSet<Conversion>();
		for (BioPAXElement ele : match.getVariables())
		{
			if (ele instanceof Conversion)
			{
				set.add((Conversion) ele);
			}
		}
		return set;
	}

	private static Integer hashSum(Set<Conversion> set)
	{
		int x = 0;
		for (Conversion inter : set)
		{
			x += inter.hashCode();
		}
		return x;
	}

	
	private static boolean isAssociated(PhysicalEntity pe, boolean input, Set<Conversion> cons)
	{
		for (Match match : Searcher.search(pe, input ? INPUT_TO_CONV_PATTERN : OUTPUT_TO_CONV_PATTERN))
		{
			if (cons.contains(match.get(1))) return true;
		}
		return false;
	}	
	
	private static boolean relatedERsDifferent(PhysicalEntity pe, Set<EntityReference> ers)
	{
		for (EntityReference er : Searcher.searchAndCollect(
			pe, RELATED_ER_PATTERN, 2, EntityReference.class))
		{
			if (!ers.contains(er)) return true;
		}
		return false;
	}
	
	// Section: Patterns
	
	private static Pattern prepareParticipantAndControllerPattern()
	{
		Pattern p = new Pattern(PhysicalEntity.class, "PE");
		p.add(ConBox.participatesInInter(), "PE", "Inter");
		p.add(ConBox.peToControl(), "PE", "Ctrl");
		p.add(ConBox.controlled(), "Ctrl", "Controlled");
		p.add(new Equality(true), "Inter", "Controlled");
		return p;
	}

	private static Pattern prepareInputToConv()
	{
		Pattern p = new Pattern(PhysicalEntity.class, "PE");
		p.add(new ParticipatesInConv(RelType.INPUT, true), "PE", "Conv");
		return p;
	}

	private static Pattern prepareOutputToConv()
	{
		Pattern p = new Pattern(PhysicalEntity.class, "PE");
		p.add(new ParticipatesInConv(RelType.OUTPUT, true), "PE", "Conv");
		return p;
	}

	private static Pattern prepareSimpleMemberPattern()
	{
		Pattern p = new Pattern(PhysicalEntity.class, "PE");
		p.add(ConBox.withSimpleMembers(), "PE", "SPE");
		p.add(new Type(SimplePhysicalEntity.class), "SPE");
		return p;
	}

	private static Pattern prepareRelatedERPattern()
	{
		Pattern p = new Pattern(PhysicalEntity.class, "PE");
		p.add(new LinkedPE(LinkedPE.Type.TO_MEMBER), "PE", "SPE");
		p.add(ConBox.peToER(), "SPE", "ER");
		return p;
	}

	private static Pattern prepareBothInputAndOutputPattern()
	{
		Pattern p = new Pattern(Conversion.class, "Conv");
		p.add(ConBox.left(), "Conv", "PE left");
		p.add(ConBox.isHuman(), "PE left");
		p.add(ConBox.right(), "PE left");
		return p;
	}

	private static Pattern prepareLoop(int loopLength)
	{
		assert loopLength > 1;
		
		Pattern p = new Pattern(ProteinReference.class, "PR1");
		p.add(ConBox.isHuman(), "PR1");
		p.add(ConBox.erToPE(), "PR1", "SPE1 in");
		p.add(new LinkedPE(LinkedPE.Type.TO_COMPLEX), "SPE1 in", "PE1 in");
		p.add(new ParticipatesInConv(RelType.INPUT, true), "PE1 in", "Conv1");
		p.add(ConBox.notControlled(), "Conv1");
		p.add(new OtherSide(), "PE1 in", "Conv", "PE1 out");
		p.add(new Equality(false), "PE1 in", "PE1 out");
		p.add(new LinkedPE(LinkedPE.Type.TO_MEMBER), "PE1 out", "SPE1 out");
		p.add(ConBox.peToER(), "SPE1 out", "PR1");

		for (int i = 2; i <= loopLength; i++)
		{
			p.add(new ParticipatesInConv(RelType.INPUT, true), "PE" + (i-1) + " out", "Conv" + i);
			p.add(ConBox.notControlled(), "Conv" + i);

			for (int j = 1; j < i; j++)
			{
				p.add(new Equality(false), "Conv" + j, "Conv" + i);
			}

			if (i == loopLength)
			{
				p.add(new OtherSide(), "PE" + (i-1) + " out", "Conv" + i, "PE1 in");
			}
			else
			{
				p.add(new OtherSide(), "PE" + (i-1) + " out", "Conv" + i, "PE" + i + " out");

				for (int j = 1; j < i; j++)
				{
					p.add(new Equality(false), "PE" + j + " out", "PE" + i + " out");
				}

				p.add(new LinkedPE(LinkedPE.Type.TO_MEMBER), "PE" + i + " out", "SPE" + i + " out");
				p.add(ConBox.peToER(), "SPE" + i + " out", "PR");
			}
		}

		return p;
	}
}
