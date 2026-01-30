// test_neurobotanica_api.js
const API_BASE = 'https://terpene-api.contessapetrini.workers.dev';

async function testHealth() {
  console.log('Testing health endpoint...');
  try {
    const response = await fetch(`${API_BASE}/api/neurobotanica/health`);
    const data = await response.json();
    console.log('Health response:', data);
    if (data.status === 'healthy' && data.engines.includes('interactions')) {
      console.log('✅ Health test passed');
      return true;
    } else {
      console.log('❌ Health test failed');
      return false;
    }
  } catch (error) {
    console.log('❌ Health test error:', error.message);
    return false;
  }
}

async function testAnalyze() {
  console.log('Testing analyze endpoint...');
  try {
    const response = await fetch(`${API_BASE}/api/neurobotanica/analyze`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        compound_ids: ['cbd', 'thc'],
        demographics: { age: 30, condition: 'weight_management' },
        customer_tier: 'computational_only'
      })
    });
    const data = await response.json();
    console.log('Analyze response:', data);
    if (data.interactions && data.bias_correction && data.synergy) {
      console.log('✅ Analyze test passed');
      return true;
    } else {
      console.log('❌ Analyze test failed');
      return false;
    }
  } catch (error) {
    console.log('❌ Analyze test error:', error.message);
    return false;
  }
}

async function runTests() {
  const healthOk = await testHealth();
  const analyzeOk = await testAnalyze();

  console.log('\nTest Results:');
  console.log('Health check:', healthOk ? '✅ PASS' : '❌ FAIL');
  console.log('Analyze endpoint:', analyzeOk ? '✅ PASS' : '❌ FAIL');

  if (healthOk && analyzeOk) {
    console.log('\n🎉 All tests passed! NeuroBotanica API integration is working.');
  } else {
    console.log('\n❌ Some tests failed. Check the API deployment.');
  }
}

runTests();