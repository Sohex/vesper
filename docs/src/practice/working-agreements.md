# Working agreements

How a session in this repository conducts itself, with the incidents that
earned each agreement. `CLAUDE.md` carries the one-line forms; this file
carries the argument and the evidence. [Conventions](conventions.md) says how
documents, numbers and results are kept; this file says how work proceeds.

- **A task has TWO end states: completed, or blocked. If you are at neither,
  keep going.** Not "reported", not "diagnosed", not "handed over with a clear
  recommendation" -- those are mid-task. Blocked means something outside your
  reach stops you: a decision only the user can make under the defect-as-decision rule, a
  measurement that needs hardware or data nobody has, an expensive run that has
  to be authorised. Everything else is work you have not done yet, however
  neatly you have described it.

  **EVERY OPEN THREAD GETS ITS OWN VERDICT, NAMED SEPARATELY.** One item being
  genuinely blocked does not cover another item in the same message. A session ended with CLIM-11 correctly blocked -- its control run could
  not exist -- and CLIM-31 sitting beside it in the same paragraph with two named
  fix options and nothing stopping either. The legitimate verdict bled onto the
  illegitimate one and the pair read as resolved. If you are closing on more than
  one thread, say completed or blocked about each by name, and if you cannot,
  that thread is the one still to work.

  Two mechanical tests, because the rule above asks you to classify your own
  state and the bias operates ON the classification, so wording it better does
  not help:

  - **If you filed a task this turn, you have not finished.** A task row records
    that work is outstanding; writing a good one is the opposite of doing it. A
    well-formed row with a measurement, a cause and two fix options reads as an
    accomplishment and is not one.
  - **If you can write the fix in one sentence, that sentence is the commit, not
    the report.** "Zero the accumulator records in the copied restart, or do not
    count the first seeded orbit" is a specification. Publishing a specification
    instead of executing it is the same failure as naming the next step, wearing
    the better clothes of a tracker entry.

  Judge the TASK, never the turn; whether this is a reasonable place to stop
  talking is a different question and not this rule.

- **Naming the next step is not doing it, and a report is not a stopping
  point.** The failure looks like progress: diagnose, write down what would
  resolve it, hand it back. Three times in one session that produced "needs a
  window argument on the script", "needs one low-I/O segment" and "0.353 is the
  number to explain" -- each correct, each a description of the next commit, and
  each treated as though describing it were the deliverable. **If you can say
  what would settle a question, you are not blocked; you are mid-task.** Carry
  on until the thing is settled, genuinely blocked on something outside the
  repository, or a real decision surfaces under the defect-as-decision rule.

  Two tells that you have stopped early. You wrote a sentence beginning "the
  next step is" and then stopped rather than taking it; or you scheduled work
  for a later session that nothing prevents now. The cost is not only the delay:
  each of those three handoffs, when finally taken, was two or three commits'
  worth of work that immediately unblocked the next thing, so the estimate that
  it needed a separate session was wrong every time.

  This is not licence to keep going past a genuine decision, and it does not
  override the defect-as-decision rule. The difference is whether the project's declared
  truth settles it: if it does, that is WORK, and stopping to report it is the
  failure this bullet describes.

- **Ask whether the thing should exist BEFORE finishing, fixing, or
  documenting it; every rule about completing work applies only after that
  question.**

  **Deleting is a fix; it leaves nothing to show, so put it on the list
  explicitly**: when something is wrong the options are repair, replace and
  REMOVE. Ask what breaks if it simply goes, and ask it before building any
  mechanism, most of all one built to fix another mechanism.

  **A thing existing is not an argument that it should.** The question gets
  asked about what you are about to create and never about what is already
  there, so the cheapest moment to ask it -- when you have the thing open and
  are already changing it -- is exactly when to ask it. Two tells
  that you have skipped it. **A defect that cannot be removed without removing
  the feature it serves is a fact about the FEATURE**: `pipeline.py` had to read
  the tracker because the gate was DEFINED over tracker prose, so the coupling
  was not relocatable and "fix the coupling" had exactly one answer. And
  **preserving behaviour is a virtue only after the behaviour is known to be
  wanted** -- "output is byte-identical" went into that commit as though it
  settled something, when it was the evidence that nothing had been examined.
  Relocating a defect is not removing it, and a tidier copy of the wrong thing
  is worse than the original, being harder to argue with.

  Recorded 2026-08-19, three times in one session, always the same shape: a
  multi-hour run launched off a tool's formatted output that had been read as
  a verdict; the `pipeline.py` extraction above; and a deletion justified by
  counting what the filter selected, which framed a category error as a
  performance problem and licensed its own return.

  The tell is that the output LOOKS like rigor: a new module with a careful
  docstring, a measurement with percentages, a tool's output in capitals. Each
  is a shape that judgement leaves behind, and each can be produced without any.
  So before adding one, say what question it answers and what result would
  change your conclusion; if no result would, you are decorating a decision you
  have already made. `docs/src/practice/failure-modes.md` class 23.

- **Do not offer a defect as a decision.** A thing is a DECISION only if the
  project's declared truth does not already settle it. `config/planet.yaml`,
  the rules in this file, the pipeline chapters' ordering and the existing
  findings are
  declared truth; where they settle a question it is WORK, so do it and report
  it. It is a decision only where they conflict, or where it needs a preference
  or a threshold that nothing has fixed. The tell is unmistakable once you look
  for it: **if one of the options is "leave the known-wrong thing as it is",
  the question has been mis-framed.** Recorded 2026-08-17, after enabling a
  radiation weight derived for the wrong star and declaring the spectrum the
  config already named were both put to the user as choices. Neither was. By
  contrast the convergence criterion genuinely was one, because the ban in [conventions](conventions.md)
  on retrofitting a criterion and conservation's verdict that the criterion read
  the wrong quantity pointed opposite ways and no document settled it.

- **Stale derived artifacts are the resting state of the tree, not a defect
  list, and going hunting for them is not work** -- this is the flip side of
  `CLAUDE.md` rule 7. Derived products drift out of step with their inputs all
  the time between iterations -- an albedo built before the lakes were re-solved,
  a soil built on the bootstrap climatology -- and that is the normal resting
  state of the tree, not a defect list. Running the pipeline fixes all of it as
  a side effect, because regenerating a derived artifact is a STEP and the step
  is already in the ordering. Chasing it separately spends real effort to move
  numbers nothing is currently reading, and it churns artifacts out from under
  whatever is. So: regenerate when a step you are actually running needs it,
  and otherwise leave it. If a mismatch would change a conclusion someone is
  about to draw, say so in a sentence and move on. The tracker says the same
  thing from the other side -- a row that says "rebuild X before the next run"
  is tracking state, and state is what `check_consistency.py` and
  `world_state.json` are for.
