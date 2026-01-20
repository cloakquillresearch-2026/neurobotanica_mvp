import asyncio
import os
os.environ['DATABASE_URL'] = 'sqlite:///./neurobotanica_dev.db'

from src.api.neurobotanica import analyze
from pydantic import BaseModel

class MockRequest(BaseModel):
    compound_ids: list[str]
    demographics: dict = {}
    customer_tier: str = "computational_only"
    plant_id: str = None

async def main():
    request = MockRequest(compound_ids=['cbd'], demographics={'age': 30}, customer_tier='computational_only', plant_id='cbd')

    try:
        result = await analyze(request)
        print('Success:', result)
    except Exception as e:
        print('Error:', e)
        import traceback
        traceback.print_exc()

asyncio.run(main())