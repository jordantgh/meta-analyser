import openai

openai.api_key = "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"

message_input = [
    {
      "role": "system", "content": "You are a helpful and intelligent assistant."
    },
    {
      "role": "user",
      "content": 
         f"""
           There are three sisters in a room. Anna is reading a book. Alice  is playing a game of chess. What is the third sister, Amanda,  doing?
         """
    }
  ]

response = openai.ChatCompletion.create(
            model="gpt-4-0613",
            messages=message_input,
            temperature=1
        )

print(response.choices[0]["message"]["content"])
