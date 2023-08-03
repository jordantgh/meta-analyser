from typing import List
from pydantic import BaseModel

class StepByStepAIResponse(BaseModel):
    title: str
    steps: List[str]
    
schema = StepByStepAIResponse.schema() # returns a dict like JSON schema

# schema content looks like below

"""
{
    'title': 'StepByStepAIResponse',
    'type': 'object',
    'properties': {'title': {'title': 'Title', 'type': 'string'},
    'steps': {'title': 'Steps', 'type': 'array', 'items': {'type': 'string'}}},
    'required': ['title', 'steps']
}
"""

import openai
import json

openai.api_key = "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"

response = openai.ChatCompletion.create(
    model="gpt-3.5-turbo-0613",
    messages=[
       {"role": "user", "content": "Explain how to assemble a PC"}
    ],
    functions=[
        {
          "name": "get_answer_for_user_query",
          "description": "Get user answer in series of steps",
          "parameters": StepByStepAIResponse.schema()
        }
    ],
    function_call={"name": "get_answer_for_user_query"}
)

output = json.loads(response.choices[0]["message"]["function_call"]["arguments"])

# output content
"""
{
    'title': 'Steps to assemble a PC',
    'steps': [
        '1. Gather all necessary components',
        '2. Prepare the PC case',
        '3. Install the power supply',
        '4. Mount the motherboard',
        '5. Install the CPU and CPU cooler',
        '6. Install RAM modules',
        '7. Install storage devices',
        '8. Install the graphics card',
        '9. Connect all cables',
        '10. Test the PC'
    ]
}
"""