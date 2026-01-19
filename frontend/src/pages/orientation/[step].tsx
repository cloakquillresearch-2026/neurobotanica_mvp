import { useRouter } from 'next/router';
import { useEffect, useState } from 'react';

export default function OrientationStep() {
  const router = useRouter();
  const { step } = router.query;
  const [answers, setAnswers] = useState<Record<string, string>>({});

  const currentStep = step === 'welcome' ? 0 : 
                      step === 'module-1' ? 1 :
                      step === 'module-2' ? 2 :
                      step === 'module-3' ? 3 :
                      step === 'complete' ? 4 : 0;

  const navigateNext = () => {
    const nextSteps = ['module-1', 'module-2', 'module-3', 'complete'];
    if (currentStep < 4) {
      router.push(`/orientation/${nextSteps[currentStep]}`);
    }
  };

  const handleQuizAnswer = (questionId: string, answer: string) => {
    setAnswers({ ...answers, [questionId]: answer });
  };

  if (!step) return null;

  return (
    <div className="min-h-screen bg-gradient-to-b from-green-50 to-white">
      <div className="max-w-4xl mx-auto px-6 py-12">
        
        {/* Progress Bar */}
        {currentStep < 4 && (
          <div className="mb-8">
            <div className="flex items-center justify-between mb-2">
              <span className="text-sm text-gray-600">Progress</span>
              <span className="text-sm text-gray-600">{currentStep} of 4</span>
            </div>
            <div className="w-full bg-gray-200 rounded-full h-2">
              <div 
                className="bg-green-600 h-2 rounded-full transition-all duration-300"
                style={{ width: `${(currentStep / 4) * 100}%` }}
              />
            </div>
          </div>
        )}

        {/* Page Content */}
        {step === 'welcome' && <WelcomePage onNext={navigateNext} />}
        {step === 'module-1' && <Module1Page onNext={navigateNext} answers={answers} onAnswer={handleQuizAnswer} />}
        {step === 'module-2' && <Module2Page onNext={navigateNext} answers={answers} onAnswer={handleQuizAnswer} />}
        {step === 'module-3' && <Module3Page onNext={navigateNext} answers={answers} onAnswer={handleQuizAnswer} />}
        {step === 'complete' && <CompletePage />}
      </div>
    </div>
  );
}

function WelcomePage({ onNext }: { onNext: () => void }) {
  return (
    <div className="text-center space-y-8">
      <h1 className="text-4xl font-bold text-gray-900">
        Welcome to NeuroBotanica Budtender Orientation
      </h1>
      <p className="text-xl text-gray-600">
        Learn how clinical evidence powers personalized cannabis recommendations
      </p>

      <div className="bg-white rounded-lg shadow-lg p-8 my-8">
        <div className="aspect-video bg-black rounded-lg overflow-hidden mb-6">
          <iframe 
            src="https://www.canva.com/design/DAG-oazPd7M/Iyp4VI6UeyLtaYNBfAbyEQ/watch?embed"
            allow="fullscreen"
            className="w-full h-full"
          />
        </div>
        
        <div className="text-left space-y-4 text-gray-700">
          <p className="text-lg">In this orientation, you will:</p>
          <ul className="list-disc list-inside space-y-2 ml-4">
            <li>Understand the NeuroBotanica mission and value proposition</li>
            <li>Learn Nevada compliance guardrails</li>
            <li>See how terpene and biomarker pairing shortens consultations</li>
            <li>Master three key talking points for customers</li>
          </ul>
          <p className="text-sm text-gray-600 pt-4">
            <strong>Estimated time:</strong> 8-10 minutes
          </p>
        </div>
      </div>

      <button
        onClick={onNext}
        className="bg-green-600 text-white px-8 py-4 rounded-lg text-lg font-semibold hover:bg-green-700 transition"
      >
        Begin Orientation
      </button>
    </div>
  );
}

function Module1Page({ onNext, answers, onAnswer }: { onNext: () => void; answers: Record<string, string>; onAnswer: (id: string, answer: string) => void }) {
  const questions = [
    { id: 'q1', question: 'What does NeuroBotanica use to power recommendations?', options: ['Personal opinions', 'Clinical evidence from 505+ studies', 'Random selection', 'Customer reviews'], correct: 1 },
    { id: 'q2', question: 'How much faster is validation with NeuroBotanica?', options: ['50%', '70%', '92%', '100%'], correct: 2 },
    { id: 'q3', question: 'What is a TS-PS-001 panel?', options: ['A type of cannabis strain', 'An evidence report showing terpene effects', 'A customer feedback form', 'A dispensary license'], correct: 1 }
  ];

  const allAnswered = questions.every(q => answers[q.id] !== undefined);
  const allCorrect = questions.every(q => answers[q.id] === q.correct.toString());
  const canContinue = allAnswered && allCorrect;

  return (
    <div className="space-y-8">
      <div className="text-center mb-8">
        <p className="text-sm text-green-600 font-semibold mb-2">Module 1 of 3</p>
        <h2 className="text-3xl font-bold text-gray-900">NeuroBotanica in 90 Seconds</h2>
      </div>

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-6">
        <h3 className="text-xl font-semibold text-gray-900">Learning Objectives</h3>
        <ul className="list-disc list-inside space-y-2 text-gray-700 ml-4">
          <li>Understand the NeuroBotanica mission and value proposition</li>
          <li>Learn Nevada compliance guardrails</li>
          <li>See how terpene and biomarker pairing shortens consultations</li>
          <li>Remember three talking points for customers</li>
        </ul>
      </div>

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-6">
        <h3 className="text-lg font-semibold text-gray-900">1. Watch the Platform Overview</h3>
        <div className="aspect-video bg-black rounded-lg overflow-hidden">
          <iframe 
            src="https://www.canva.com/design/DAG-1vK9JMQ/wHgVgLi8sMW2-Ldy9_FMlQ/watch?embed"
            allow="fullscreen"
            className="w-full h-full"
          />
        </div>
      </div>

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-6">
        <h3 className="text-lg font-semibold text-gray-900">2. Review Key Statistics</h3>
        <div className="grid grid-cols-3 gap-4">
          <div className="bg-green-50 p-4 rounded-lg text-center">
            <p className="text-3xl font-bold text-green-600">505+</p>
            <p className="text-sm text-gray-600">Clinical Studies</p>
          </div>
          <div className="bg-green-50 p-4 rounded-lg text-center">
            <p className="text-3xl font-bold text-green-600">92%</p>
            <p className="text-sm text-gray-600">Faster Validation</p>
          </div>
          <div className="bg-green-50 p-4 rounded-lg text-center">
            <p className="text-3xl font-bold text-green-600">85%+</p>
            <p className="text-sm text-gray-600">Accuracy</p>
          </div>
        </div>
      </div>

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-6">
        <h3 className="text-lg font-semibold text-gray-900">3. Commitment Check</h3>
        <p className="text-gray-700">Can you explain the TS-PS-001 value to a customer in one sentence?</p>
        <textarea
          className="w-full p-4 border border-gray-300 rounded-lg"
          rows={3}
          placeholder="Type your answer here..."
          value={answers['commitment'] || ''}
          onChange={(e) => onAnswer('commitment', e.target.value)}
        />
      </div>

      <QuizSection
        title="Quick Check"
        questions={questions}
        answers={answers}
        onAnswer={onAnswer}
      />

      {!canContinue && allAnswered && (
        <div className="bg-red-50 border-l-4 border-red-600 p-4">
          <p className="text-red-800">Please answer all questions correctly before continuing.</p>
        </div>
      )}

      <button
        onClick={onNext}
        disabled={!canContinue}
        className={`w-full px-8 py-4 rounded-lg text-lg font-semibold transition ${
          canContinue 
            ? 'bg-green-600 text-white hover:bg-green-700 cursor-pointer' 
            : 'bg-gray-300 text-gray-500 cursor-not-allowed'
        }`}
      >
        Continue to Module 2
      </button>
    </div>
  );
}

function Module2Page({ onNext, answers, onAnswer }: { onNext: () => void; answers: Record<string, string>; onAnswer: (id: string, answer: string) => void }) {
  const personas = [
    { name: 'Sarah, 45', condition: 'Chronic Pain', focus: 'β-caryophyllene (non-psychoactive relief)', method: 'Topicals & tinctures', concerns: 'Wants to avoid THC, needs daytime relief' },
    { name: 'Marcus, 32', condition: 'Anxiety', focus: 'Linalool (calming effects)', method: 'Edibles preferred', concerns: 'Worried about paranoia, needs evening relaxation' },
    { name: 'Janet, 58', condition: 'Nausea from Chemotherapy', focus: 'Limonene + ginger adjuvant', method: 'Sublingual tinctures', concerns: 'Sensitive stomach, needs fast relief' }
  ];

  const questions = [
    { id: 'q4', question: 'Which terpene would you recommend for Sarah\'s chronic pain?', options: ['Linalool', 'β-caryophyllene', 'Myrcene', 'Limonene'], correct: 1 },
    { id: 'q5', question: 'What\'s the key concern for Marcus with anxiety?', options: ['Nausea', 'Paranoia', 'Pain', 'Insomnia'], correct: 1 },
    { id: 'q6', question: 'Which adjuvant helps with Janet\'s nausea?', options: ['Black pepper', 'MCT oil', 'Ginger', 'Turmeric'], correct: 2 }
  ];

  const allAnswered = questions.every(q => answers[q.id] !== undefined);
  const allCorrect = questions.every(q => answers[q.id] === q.correct.toString());
  const canContinue = allAnswered && allCorrect;

  return (
    <div className="space-y-8">
      <div className="text-center mb-8">
        <p className="text-sm text-green-600 font-semibold mb-2">Module 2 of 3</p>
        <h2 className="text-3xl font-bold text-gray-900">Case Study Sprint</h2>
      </div>

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-6">
        <h3 className="text-xl font-semibold text-gray-900">Learning Objectives</h3>
        <ul className="list-disc list-inside space-y-2 text-gray-700 ml-4">
          <li>Map customer conditions to terpene/adjuvant narratives</li>
          <li>Practice phrasing recommendations with clinical language</li>
          <li>Capture confidence ratings for analytics</li>
        </ul>
      </div>

      <div className="space-y-4">
        <h3 className="text-xl font-semibold text-gray-900">Interactive Personas</h3>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          {personas.map((persona, idx) => (
            <div key={idx} className="bg-white rounded-lg shadow-lg p-6 border-2 border-gray-200 hover:border-green-500 cursor-pointer transition">
              <h4 className="text-lg font-bold text-gray-900 mb-2">{persona.name}</h4>
              <p className="text-sm text-gray-600 mb-1"><strong>Condition:</strong> {persona.condition}</p>
              <p className="text-sm text-gray-600 mb-1"><strong>Focus:</strong> {persona.focus}</p>
              <p className="text-sm text-gray-600 mb-1"><strong>Method:</strong> {persona.method}</p>
              <p className="text-sm text-gray-600"><strong>Concerns:</strong> {persona.concerns}</p>
            </div>
          ))}
        </div>
      </div>

      <QuizSection
        title="Scenario Questions"
        questions={questions}
        answers={answers}
        onAnswer={onAnswer}
      />

      {!canContinue && allAnswered && (
        <div className="bg-red-50 border-l-4 border-red-600 p-4">
          <p className="text-red-800">Please answer all questions correctly before continuing.</p>
        </div>
      )}

      <button
        onClick={onNext}
        disabled={!canContinue}
        className={`w-full px-8 py-4 rounded-lg text-lg font-semibold transition ${
          canContinue 
            ? 'bg-green-600 text-white hover:bg-green-700 cursor-pointer' 
            : 'bg-gray-300 text-gray-500 cursor-not-allowed'
        }`}
      >
        Continue to Module 3
      </button>
    </div>
  );
}

function Module3Page({ onNext, answers, onAnswer }: { onNext: () => void; answers: Record<string, string>; onAnswer: (id: string, answer: string) => void }) {
  const questions = [
    { id: 'q7', question: 'What is the entourage effect?', options: ['A type of terpene', 'Synergy between cannabinoids and terpenes', 'A cannabis strain', 'An extraction method'], correct: 1 },
    { id: 'q8', question: 'Which terpene is best for anxiety?', options: ['Myrcene', 'Linalool', 'Pinene', 'Humulene'], correct: 1 },
    { id: 'q9', question: 'What does piperine (black pepper) do?', options: ['Adds flavor', 'Increases cannabinoid absorption', 'Reduces THC effects', 'Nothing'], correct: 1 },
    { id: 'q10', question: 'What is β-caryophyllene unique for?', options: ['It smells like lemons', 'It acts on CB2 receptors', 'It causes sedation', 'It increases appetite'], correct: 1 },
    { id: 'q11', question: 'Where do you find terpene data?', options: ['On the product label', 'Certificate of Analysis (CoA)', 'Customer reviews', 'Budtender opinion'], correct: 1 }
  ];

  const allAnswered = questions.every(q => answers[q.id] !== undefined);
  const allCorrect = questions.every(q => answers[q.id] === q.correct.toString());
  const canContinue = allAnswered && allCorrect;

  return (
    <div className="space-y-8">
      <div className="text-center mb-8">
        <p className="text-sm text-green-600 font-semibold mb-2">Module 3 of 3</p>
        <h2 className="text-3xl font-bold text-gray-900">Terpenes & Adjuvants</h2>
      </div>

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-6">
        <h3 className="text-xl font-semibold text-gray-900">Learning Objectives</h3>
        <ul className="list-disc list-inside space-y-2 text-gray-700 ml-4">
          <li>Understand what terpenes are and their therapeutic effects</li>
          <li>Learn the entourage effect and why whole-plant extracts work better</li>
          <li>Know key terpenes: linalool, β-caryophyllene, myrcene, limonene</li>
          <li>Understand adjuvants and their role in enhancing efficacy</li>
          <li>Apply knowledge to read a terpene profile</li>
        </ul>
      </div>

      <ContentSection
        number={1}
        title="What are Terpenes?"
        content="Terpenes are aromatic compounds found in plants that give them their distinctive scents and flavors. In cannabis, terpenes work synergistically with cannabinoids to produce therapeutic effects."
        tip="Think of terpenes like essential oils—they give plants their scent and medicinal qualities."
      />

      <ContentSection
        number={2}
        title="The Entourage Effect"
        content="The entourage effect describes how cannabinoids and terpenes work together synergistically to enhance therapeutic benefits. For example, myrcene increases THC absorption across the blood-brain barrier."
        tip="Whole-plant extracts often work better than isolated compounds because of this synergy."
      />

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-4">
        <h3 className="text-lg font-semibold text-gray-900">3. Key Terpenes to Know</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <TerpeneCard name="Linalool" effects="Calming, anti-anxiety" />
          <TerpeneCard name="β-caryophyllene" effects="Anti-inflammatory, CB2 agonist" />
          <TerpeneCard name="Myrcene" effects="Sedation, muscle relaxation" />
          <TerpeneCard name="Limonene" effects="Mood elevation, anti-anxiety" />
        </div>
      </div>

      <ContentSection
        number={4}
        title="What are Adjuvants?"
        content="Adjuvants are substances that enhance the absorption or effectiveness of cannabinoids and terpenes. Examples include piperine (black pepper) which increases cannabinoid absorption, MCT oil for fat-soluble delivery, and ginger for anti-nausea synergy."
      />

      <ContentSection
        number={5}
        title="Practical Application"
        content="Learn to read lab terpene profiles by identifying dominant terpenes (usually listed in mg/g or %), matching them to customer goals, and explaining benefits in simple terms. Always reference the lab certificate of analysis (CoA)."
      />

      <QuizSection
        title="Final Assessment"
        questions={questions}
        answers={answers}
        onAnswer={onAnswer}
      />

      {!canContinue && allAnswered && (
        <div className="bg-red-50 border-l-4 border-red-600 p-4">
          <p className="text-red-800">Please answer all questions correctly before continuing.</p>
        </div>
      )}

      <button
        onClick={onNext}
        disabled={!canContinue}
        className={`w-full px-8 py-4 rounded-lg text-lg font-semibold transition ${
          canContinue 
            ? 'bg-green-600 text-white hover:bg-green-700 cursor-pointer' 
            : 'bg-gray-300 text-gray-500 cursor-not-allowed'
        }`}
      >
        Complete Orientation
      </button>
    </div>
  );
}

function CompletePage() {
  const router = useRouter();

  const handleComplete = () => {
    // Mark orientation as complete in localStorage
    if (typeof window !== 'undefined') {
      localStorage.setItem('orientation_complete', 'true');
    }
    router.push('/');
  };

  return (
    <div className="text-center space-y-8">
      <div className="text-6xl mb-4">🎉</div>
      <h1 className="text-4xl font-bold text-gray-900">
        Congratulations! Orientation Complete
      </h1>
      <p className="text-xl text-gray-600">
        You've completed all three modules and are ready to use NeuroBotanica at the counter.
      </p>

      <div className="bg-white rounded-lg shadow-lg p-8 space-y-4">
        <div className="flex items-center justify-center space-x-4 mb-6">
          <div className="flex items-center space-x-2">
            <div className="w-8 h-8 bg-green-600 rounded-full flex items-center justify-center text-white">✓</div>
            <span className="text-gray-700">Module 1</span>
          </div>
          <div className="flex items-center space-x-2">
            <div className="w-8 h-8 bg-green-600 rounded-full flex items-center justify-center text-white">✓</div>
            <span className="text-gray-700">Module 2</span>
          </div>
          <div className="flex items-center space-x-2">
            <div className="w-8 h-8 bg-green-600 rounded-full flex items-center justify-center text-white">✓</div>
            <span className="text-gray-700">Module 3</span>
          </div>
        </div>

        <div className="text-left space-y-3 text-gray-700">
          <h3 className="font-semibold text-lg">Next Steps:</h3>
          <ul className="list-disc list-inside space-y-2 ml-4">
            <li>Your progress syncs automatically to the Budtender tablet</li>
            <li>You now have access to Sandbox Practice mode and Certification Quizzes</li>
            <li>Start using NeuroBotanica with real customers</li>
          </ul>
        </div>
      </div>

      <div className="flex flex-col sm:flex-row gap-4 justify-center">
        <button
          onClick={() => router.push('/?sandbox=1')}
          className="bg-green-600 text-white px-8 py-4 rounded-lg text-lg font-semibold hover:bg-green-700 transition"
        >
          Launch Sandbox Practice
        </button>
        <button
          onClick={handleComplete}
          className="bg-white border-2 border-green-600 text-green-600 px-8 py-4 rounded-lg text-lg font-semibold hover:bg-green-50 transition"
        >
          Return to Budtender App
        </button>
      </div>
    </div>
  );
}

// Helper Components
function ContentSection({ number, title, content, tip }: { number: number; title: string; content: string; tip?: string }) {
  return (
    <div className="bg-white rounded-lg shadow-lg p-8 space-y-4">
      <h3 className="text-lg font-semibold text-gray-900">{number}. {title}</h3>
      <p className="text-gray-700">{content}</p>
      {tip && (
        <div className="bg-green-50 border-l-4 border-green-600 p-4">
          <p className="text-sm text-green-800"><strong>💡 Tip:</strong> {tip}</p>
        </div>
      )}
    </div>
  );
}

function TerpeneCard({ name, effects }: { name: string; effects: string }) {
  return (
    <div className="bg-green-50 p-4 rounded-lg">
      <h4 className="font-semibold text-gray-900 mb-1">{name}</h4>
      <p className="text-sm text-gray-600">{effects}</p>
    </div>
  );
}

function QuizSection({ title, questions, answers, onAnswer }: { 
  title: string; 
  questions: Array<{ id: string; question: string; options: string[]; correct: number }>; 
  answers: Record<string, string>; 
  onAnswer: (id: string, answer: string) => void 
}) {
  return (
    <div className="bg-white rounded-lg shadow-lg p-8 space-y-6">
      <h3 className="text-lg font-semibold text-gray-900">{title}</h3>
      {questions.map((q, idx) => (
        <div key={q.id} className="space-y-3">
          <p className="font-medium text-gray-900">{idx + 1}. {q.question}</p>
          <div className="space-y-2">
            {q.options.map((option, optIdx) => (
              <label key={optIdx} className="flex items-center space-x-3 cursor-pointer">
                <input
                  type="radio"
                  name={q.id}
                  value={optIdx.toString()}
                  checked={answers[q.id] === optIdx.toString()}
                  onChange={(e) => onAnswer(q.id, e.target.value)}
                  className="w-4 h-4 text-green-600"
                />
                <span className="text-gray-700">{option}</span>
              </label>
            ))}
          </div>
          {answers[q.id] && (
            <p className={`text-sm font-semibold ${answers[q.id] === q.correct.toString() ? 'text-green-600' : 'text-red-600'}`}>
              {answers[q.id] === q.correct.toString() ? '✓ Correct!' : '✗ Incorrect. Please try again.'}
            </p>
          )}
        </div>
      ))}
    </div>
  );
}