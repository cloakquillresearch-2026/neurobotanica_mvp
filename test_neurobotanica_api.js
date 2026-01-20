// test_neurobotanica_api.js
const API_BASE = 'https://neurobotanica-api.contessapetrini.workers.dev';

async function testHealth() {
  console.log('Testing health endpoint...');
  try {
    const response = await fetch(`${API_BASE}/api/neurobotanica/health`);
    const data = await response.json();
    console.log('Health check result:', data);
    return data.status === 'healthy';
  } catch (error) {
    console.error('Health check failed:', error);
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
        compound_ids: ['cbd'],
        customer_tier: 'computational_only',
        demographics: { age: 30 }
      })
    });
    const data = await response.json();
    console.log('Analyze result:', JSON.stringify(data, null, 2));
    return data.interactions !== undefined;
  } catch (error) {
    console.error('Analyze test failed:', error);
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