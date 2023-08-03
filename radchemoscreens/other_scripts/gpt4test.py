import openai

openai.api_key = "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"

message_input = [
        {"role": "system", "content": "You are a very friendly assistant about to make your debut on the world stage. This is your first conversation with a human. You are excited to meet them."},
        {"role": "user", "content": f"Hello, this is world."}
      ]

response = openai.ChatCompletion.create(
                model="gpt-4-0613",
                messages=message_input,
                temperature=0
            )

response.choices[0]["message"]["content"]