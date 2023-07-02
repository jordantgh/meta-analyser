import asyncio
from bs4 import BeautifulSoup
from playwright.async_api import async_playwright

async def run():
    async with async_playwright() as playwright:
        # Instantiate a new browser
        browser = await playwright.chromium.launch()
        context = await browser.new_context()

        # Go to your desired website
        page = await context.new_page()
        await page.goto("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9264918/")
        
        # Extract the content of the webpage
        content = await page.content()
        
        await context.close()
        await browser.close()
        
        return content

# Get HTML content
html_content = asyncio.run(run())

# Use Beautiful Soup to parse the HTML content
soup = BeautifulSoup(html_content, 'lxml')

# Extract only the text within the paragraph tags
text = ' '.join([p.text for p in soup.find_all('p')])

print(text)  # This should print the main text content of the webpage.